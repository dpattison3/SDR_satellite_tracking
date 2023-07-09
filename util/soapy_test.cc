#include <SoapySDR/Version.hpp>
#include <SoapySDR/Modules.hpp>
#include <SoapySDR/Registry.hpp>
#include <SoapySDR/Device.hpp>
#include <SoapySDR/ConverterRegistry.hpp>

#include <AudioFile.h>

#include <iostream>
#include <complex>
#include <chrono>

#include "util.h"


// It's assumed (and verified) that there is only one receiving channel.
constexpr size_t RX_CHANNEL_INDEX = 0;


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
        print(setting_info.description, setting_info.value, setting_value);
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
    const double current_sample_rate = device->getSampleRate(SOAPY_SDR_RX, 0);
    print("\ncurrent sample rate (Mhz):", current_sample_rate / 1e6);
}


int main() {
    // NOTE(dominic): I've written this mostly to interact with the SDRPlay RSPdx. A couple general
    // notes on that specifically that I've found and don't want to forget.
    // The gain setting is confusingly accessible in both the
    // settings API (getSettingInfo, readSetting, writeSetting) and the gain API (getGain, setGain).
    // The settings interface is confusing because it lists a set of possible options that do not
    // match the gain API. Some of the options don't seem to make any sense, but the options don't
    // appear to be checked at all. I've confirmed it's updating the same low-level value. Prefer
    // using the gain API directly instead.
    //
    // Only Antenna B supports bias T capabilities. This isn't checked by the SW, so you could
    // "enable" but it will have no affect.

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

    // TODO: provide a way to change the current settings.
    PrintSettingsInfo(device);

    // Prep for reading data from the antenna.
    device->setFrequency(SOAPY_SDR_RX, RX_CHANNEL_INDEX, 103.7e6);

    // format: "CF32" -  complex float32 (8 bytes per element)
    SoapySDR::Stream* const stream = device->setupStream(SOAPY_SDR_RX, /* format */ "CF32");
    const size_t max_buffer_size = device->getStreamMTU(stream);

    // Allocate buffer. Kinda gross memory management. See if we can update to a unique_ptr under
    // the hood.
    void *buffs[1] = { nullptr };
    buffs[0] = malloc(max_buffer_size * 2 * sizeof(float));

    // output data
    // std::vector<std::complex<float>> raw_output_data;
    // raw_output_data.reserve(max_buffer_size);
    std::complex<float> current_raw_iq, previous_raw_iq;
    float demodulated_value;

    std::vector<float> audio_output;
    audio_output.reserve(80000);

    std::cout << "max buffer size: " << max_buffer_size << std::endl;
    device->activateStream(stream);




    int num_iterations = 500;
    std::vector<int64_t> timing;
    std::vector<double> actual_sample_rates;
    timing.reserve(num_iterations);
    actual_sample_rates.reserve(num_iterations);

    long long read_time_ns = 0, previous_read_time_ns = 0;
    int read_stream_return_flags = 0;
    int total_num_read = 0;

std::chrono::steady_clock::time_point timer, previous_timer;
std::chrono::steady_clock::time_point end;
std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for (int outer_loop = 0; outer_loop < num_iterations; outer_loop++) {
        const int num_read = device->readStream(stream, buffs, max_buffer_size,
            read_stream_return_flags, read_time_ns);
        total_num_read += num_read;

        // std::cout << "Read: " << num_read << std::endl;

        // Copy the data out of the raw buffer.
        // Code copied from elsewhere. The void** is an array of complex floating point numbers.
        // The array is of form:
        // [r0, i0, r1, i1], where r{n} is the real part of the n th sample and c{n} is the imaginary
        // part of the n th sample.
        float *stream_output = (float *)buffs[0];

        // std::cout << "data" << std::endl;
        // raw_output_data.resize(num_read);
        previous_raw_iq = {stream_output[0], stream_output[1]};
        for (int i = 1; i < num_read; ++i) {
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
    end = std::chrono::steady_clock::now();

    int64_t elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    double average_sample_rate = 1e9 * total_num_read / elapsed_time;
std::cout << "Time difference = " << elapsed_time << "[µs]" << std::endl << "avg sample rate: " <<  average_sample_rate << std::endl;


    device->closeStream(stream);
    free(buffs[0]);


    // std::cout << std::endl << std::endl << "sample rates:" << std::endl;
    // for (const double rate : actual_sample_rates) {
    //     std::cout << rate << std::endl;
    // }



    AudioFile<float> wav_file;
    wav_file.setNumChannels(1);
    const int audio_sample_rate = 44100;
    wav_file.setNumSamplesPerChannel(audio_output.size());
    for (int i = 0; i < wav_file.getNumSamplesPerChannel(); ++i) {
        wav_file.samples[/* channel */ 0][i] = audio_output[i];
    }
    wav_file.save("data/out.wav", AudioFileFormat::Wave);

    print("audio length:", audio_output.size(), " data size:", sizeof(float) * audio_output.size() / 1e6);



    // Separately, try out the wav file writer.
    // Unfortunately, this requires pre-allocating the entire data buffer. There's no incremental
    // writing functionality. That seems pretty limiting...
    // AudioFile<float> wav_file;
    // wav_file.setNumChannels(1);
    // const int audio_sample_rate = 44100;
    // wav_file.setNumSamplesPerChannel(audio_sample_rate * 5);
    // for (int i = 0; i < wav_file.getNumSamplesPerChannel(); ++i) {
    //     wav_file.samples[/* channel */ 0][i] = std::cos(
    //         static_cast<float>(i) / static_cast<float>(audio_sample_rate) * 400 * 2 * static_cast<float>(M_PI)
    //     );
    // }
    // wav_file.save("sine.wav", AudioFileFormat::Wave);

/*
std::vector<std::string> SoapySDR::Device::listAntennas(const int, const size_t) const
{
    return std::vector<std::string>();
}

void SoapySDR::Device::setAntenna(const int, const size_t, const std::string &)
{
    return;
}

std::string SoapySDR::Device::getAntenna(const int, const size_t) const
{
    return "";
}

*/

}