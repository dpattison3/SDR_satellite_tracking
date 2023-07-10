#include <iostream>
#include <complex>
#include <chrono>

#include <AudioFile.h>
#include <SoapySDR/Version.hpp>
#include <SoapySDR/Modules.hpp>
#include <SoapySDR/Registry.hpp>
#include <SoapySDR/Device.hpp>
#include <SoapySDR/ConverterRegistry.hpp>

#include "util.h"


// Convenience structs for storing the desired settings.
// Setting allows for storing the key along with the value. These are done as strings, which is
// bad design in the API (writeSetting). It'd be better to template so that type checking is
// explicit.
struct Setting {
    std::string key;
    std::string value;
};

struct Configuration {
    // Key value pairs that are updated with the SoapySDR::Device::writeSetting API.
    std::vector<Setting> device_settings;
    // Center frequency that the SDR is tuned to.
    double tuned_frequency;           // [hz]
    // Sample rate of the SDR.
    double sample_frequency;          // [hz]
    // Sample rate of the output audio file after demodulation.
    uint32_t audio_sample_frequency;  // [hz]
    // How long to record for. Be sure this is a reasonable length that can fit in memory.
    double recording_length;          // [s]
    // Where to save the resulting wav file.
    std::string output_audio_file;
    // Where to save the resulting raw IQ file.
    std::string output_iq_file;
};

// It's assumed (and verified) that there is only one receiving channel.
constexpr size_t RX_CHANNEL_INDEX = 0;

const Configuration DEFAULT_CONFIGURATION = {
    .device_settings = {
        // IQ correction
        {"iqcorr_ctrl", "true"},
        // Bias T enabling. Note, this is only available on antenna B.
        {"biasT_ctrl", "false"},
        // RF notch filter to block AM/FM frequencies.
        {"rfnotch_ctrl", "false"},
        // DAB (digital audio broadcast) notch filter.
        {"dabnotch_ctrl", "false"},
        // AGC (automatric gain control) setpoint measured in dBfs.
        {"agc_setpoint", "-30"},
    },
    .tuned_frequency = 103.7e6,
    .sample_frequency = 2e6,
    .audio_sample_frequency = 44100,
    .recording_length = 5,
    .output_audio_file = "data/out_alt.wav",
    .output_iq_file = "data/out_iq.rawiq"
};


void PrintSoapyInfo() {
    print("SoapySDR software info:");
    print("Lib Version:", SoapySDR::getLibVersion());
    print("API Version:", SoapySDR::getAPIVersion());
    print("ABI Version:", SoapySDR::getABIVersion());
    print("Install root:", SoapySDR::getRootPath());

    SoapySDR::KwargsList devices = SoapySDR::Device::enumerate();
    print("\nDetected ", devices.size(), " devices");
    for (const auto& device : devices) {
        for (const auto& item : device) {
            print(item.first, ":", item.second);
        }
    }
}

void PrintSdrDeviceInfo(const SoapySDR::Device* const device) {
    const SoapySDR::Kwargs hw_info = device->getHardwareInfo();
    print("\nHW device info:");
    for (const auto& item : hw_info) {
        print(item.first, "->", item.second);
    }

}

void PrintAntennaInfo(const SoapySDR::Device* const device) {
    const std::vector<std::string> antenna_list = device->listAntennas(SOAPY_SDR_RX, RX_CHANNEL_INDEX);
    print("\nFound ", antenna_list.size(), " antennas available.");
    for (const std::string& antenna : antenna_list) {
        print(antenna);
    }

    print("Current antenna: ", device->getAntenna(SOAPY_SDR_RX, RX_CHANNEL_INDEX));
}


void PrintSettingsInfo(const SoapySDR::Device* const device) {
    SoapySDR::ArgInfoList settings_info = device->getSettingInfo();
    print("\nAvailable settings:");
    for (const SoapySDR::ArgInfo& setting_info : settings_info) {
        const std::string setting_value = device->readSetting(setting_info.key);
        print(setting_info.description, setting_value);
    }

    // Gain is done separately.
    const SoapySDR::Range gain_range = device->getGainRange(SOAPY_SDR_RX, RX_CHANNEL_INDEX);
    const double current_gain = device->getGain(SOAPY_SDR_RX, RX_CHANNEL_INDEX);
    const bool agc_enabled = device->getGainMode(SOAPY_SDR_RX, RX_CHANNEL_INDEX);
    print(
        "\nGain settings:\nrange:", gain_range.minimum(), "-", gain_range.maximum(), "\ngain:", current_gain,
        "\nAGC enabled:", agc_enabled
    );

    // Sample rate is done separately.
    const double current_sample_rate = device->getSampleRate(SOAPY_SDR_RX, RX_CHANNEL_INDEX);
    print("\ncurrent sample rate (Mhz):", current_sample_rate / 1e6);
}

void SetSettings(SoapySDR::Device* const device, const Configuration& configuration) {
    for (const Setting& setting : configuration.device_settings) {
        device->writeSetting(setting.key, setting.value);
    }
    
    device->setFrequency(SOAPY_SDR_RX, RX_CHANNEL_INDEX, configuration.tuned_frequency);
    device->setSampleRate(SOAPY_SDR_RX, RX_CHANNEL_INDEX, configuration.sample_frequency);
}

void RecordFMAudio(SoapySDR::Device* const device, const Configuration& configuration) {
    // format: "CF32" -  complex float32 (8 bytes per timestep)
    SoapySDR::Stream* const stream = device->setupStream(SOAPY_SDR_RX, /* format */ "CF32");
    const size_t max_buffer_size = device->getStreamMTU(stream);

    // Allocate buffer. Kinda gross memory management. See if we can update to a unique_ptr under
    // the hood.
    void *buffs[1] = { nullptr };
    buffs[0] = malloc(max_buffer_size * 2 * sizeof(float));

    print("Starting recording for ", configuration.recording_length, "seconds with a buffer size of",
          max_buffer_size);
    device->activateStream(stream);

    // We read the device in chunks into the buffer, then process that.
    const size_t num_iterations = static_cast<size_t>(
        std::ceil(configuration.sample_frequency * configuration.recording_length / max_buffer_size)
    );

    // Annoyingly, the wav file writer used here does not have the capability to write to the file
    // incrementally. Really, this is likely fine because the audio file is small enough to easily
    // fit into memory. Also, writing to disk might be slow and affect our abillity to consistently
    // stream data from device.
    std::vector<float> audio_output;
    audio_output.reserve(static_cast<size_t>(std::ceil(
        configuration.audio_sample_frequency * configuration.recording_length
    )));

    // Incrementally write to a binary file for the raw IQ data.
    auto raw_iq_file = std::fstream(configuration.output_iq_file, std::ios::out | std::ios::binary);

    // Preallocate loop varialbes.
    long long read_time_ns = 0, previous_read_time_ns = 0;
    int read_stream_return_flags = 0, num_samples_read = 0;
    std::complex<float> current_raw_iq, previous_raw_iq;
    float demodulated_value;

    // Record the start and end times as well as total samples collected. This can indicate if our
    // reading is too slow.
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    int total_num_samples_read = 0;

    print("running for", num_iterations, "iterations");
    for (int j = 0; j < num_iterations; j++) {
        num_samples_read = device->readStream(stream, buffs, max_buffer_size,
            read_stream_return_flags, read_time_ns);
        total_num_samples_read += num_samples_read;

        if (num_samples_read <= 0) {
            print_error("Invalid read from device! num_samples_read", num_samples_read);
            break;
        }

        raw_iq_file.write((char*)buffs[0], sizeof(float)*2*num_samples_read);
        raw_iq_file.seekp(0, std::ios_base::end);

        // Copy the data out of the raw buffer.
        // Code copied from elsewhere. The void** is an array of complex floating point numbers.
        // The array is of form:
        // [r0, i0, r1, i1], where r{n} is the real part of the n th sample and c{n} is the imaginary
        // part of the n th sample.
        float *stream_output = (float *)buffs[0];

        if (j == 0) {
            previous_raw_iq = {stream_output[0], stream_output[1]};
        }

        for (int i = 0; i < num_samples_read; ++i) {
            // std:: cout << stream_output[2*i] << " " << stream_output[2*i+1]  << std::endl;
            // raw_output_data[i] = {stream_output[2*i], stream_output[2*i+1]};
            current_raw_iq = {stream_output[2*i], stream_output[2*i+1]};

            // TODO: filter with liquid_firdes_kaiser
            // see src/modules/modem/analog/ModemFMStereo.cpp
            // https://liquidsdr.org/doc/firfilt/
            // https://liquidsdr.org/doc/firdes/

            demodulated_value = std::arg( std::conj(previous_raw_iq) * current_raw_iq );

            // dumb decimation...
            // do with https://liquidsdr.org/doc/msresamp/ msresamp_crcf_create
            if (i % 50 == 0) {
                audio_output.push_back(demodulated_value);
            }
            previous_raw_iq = current_raw_iq;
        }
    }

    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const double elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    const double average_sample_rate = total_num_samples_read / elapsed_time;
    print(
        "Finished reading data. Total time elapsed",
        elapsed_time,
        "seconds with an average sample rate of",
        average_sample_rate
    );

    device->closeStream(stream);
    free(buffs[0]);
    raw_iq_file.close();

    print("Saving audio file with:", audio_output.size(), "samples.");
    AudioFile<float> wav_file;
    wav_file.setNumChannels(1);
    const int audio_sample_rate = configuration.audio_sample_frequency;
    wav_file.setNumSamplesPerChannel(audio_output.size());
    for (int i = 0; i < wav_file.getNumSamplesPerChannel(); ++i) {
        wav_file.samples[/* channel */ 0][i] = audio_output[i];
    }
    wav_file.save(configuration.output_audio_file, AudioFileFormat::Wave);
}


int main() {
    // NOTE(dominic): I've written this mostly to interact with the SDRPlay RSPdx. A couple general
    // notes on that specifically that I've found and don't want to forget.
    //
    // The gain setting is confusingly accessible in both the
    // settings API (getSettingInfo, readSetting, writeSetting) and the gain API (getGain, setGain).
    // The settings interface is confusing because it lists a set of possible options that do not
    // match the gain API. Some of the options don't seem to make any sense, but the options don't
    // appear to be checked at all. I've confirmed it's updating the same low-level value. Prefer
    // using the gain API directly instead.
    //
    // Only Antenna B supports bias T capabilities. This isn't checked by the SW, so you could
    // "enable" but it will have no affect.

    // TODO: Make these settings configurable as user input.
    const Configuration configuration = DEFAULT_CONFIGURATION;

    PrintSoapyInfo();

    // Assumes there is only a single device attached. Otherwise, this will need to be provided
    // arguments.
    print("\nTrying to load device");
    SoapySDR::Device* device = SoapySDR::Device::make();
    if (device == nullptr) {
        print_error("Failed to create SDR device.");
        return 1;
    }

    print("Successfully loaded device");
    PrintSdrDeviceInfo(device);

    // There should only be a single channel.
    const size_t device_rx_channels = device->getNumChannels(SOAPY_SDR_RX);
    if (device_rx_channels != 1) {
        print_error("Found unexpected number of RX channels: ", device_rx_channels);
        return 1;
    }

    // TODO: provide a way to change the current antenna.
    PrintAntennaInfo(device);

    // TODO: provide a way to change the current settings through user input.
    SetSettings(device, configuration);
    PrintSettingsInfo(device);
    RecordFMAudio(device, configuration);
}