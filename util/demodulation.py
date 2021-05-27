import numpy as np

from util.filtering import filter_complex_signal, filter_real_signal


# Cutoff used in filtering the audio signal. Since this is the top of the audio spectrum, we just
# define it as a constant.
AUDIO_LPF_FREQUENCY = 20e3


def demodulate_signal(signal, sample_rate, band_filter_cutoff=80000, downsample_rate=25):
    """
    Demodulates an FM signal. First, perform a LPF of input to the appropriate bandwidth. Next,
    diff the phase angles to get the data signal. LPF that signal to just the audible range.
    Finally, downsample to a more reasonable sample rate.

    Args:
        signal: complex baseband input signal
        sample_rate: sample rate of raw signal
        band_filter_cutoff: the first LPF of raw signal
        downsample_rate: the output will have sampling rate of sample_rate/downsample_rate

    Returns:
        tuple of output audio signal and its associated sample rate
    """

    filtered_signal = filter_complex_signal(signal=signal, cutoff_frequency=band_filter_cutoff,
                                            sample_rate=sample_rate)

    phase_angles = np.angle(filtered_signal)
    # Difference the phase angle and convert to +/-pi range.
    phase_angle_diff = np.unwrap(np.diff(phase_angles))

    filtered_phase_angle_diff = filter_real_signal(
        signal=phase_angle_diff, cutoff_frequency=AUDIO_LPF_FREQUENCY, sample_rate=sample_rate)


    phase_angle_diff_downsampled = filtered_phase_angle_diff[::downsample_rate]
    phase_angle_diff_downsampled_sample_rate = sample_rate / float(downsample_rate)

    return phase_angle_diff_downsampled, phase_angle_diff_downsampled_sample_rate


def chunked_demodulate_signal(signal, sample_rate, chunk_size=1, band_filter_cutoff=80000,
                              downsample_rate=25):
    """
    Demodulates an FM signal. see demodulate_signal above.
    This function processes the signal in chunks to be lighter on memory usage. This will take a hit
    on speed.
    TODO(dominic): Think a bit more on if there is an off by one error in how we crop and process
    the signals. This is likely to be pretty insignificant given that it's at a high sample rate.
    There are also likely lurking issues with if any chunk is too small.

    Args:
        signal: complex baseband input signal
        sample_rate: sample rate of raw signal
        chunk_size: in seconds the amount of data that is processed at one time.
        band_filter_cutoff: the first LPF of raw signal
        downsample_rate: the output will have sampling rate of sample_rate/downsample_rate

    Returns:
        tuple of output audio signal and its associated sample rate
    """

    signal_length = len(signal)
    start_index = 0
    chunk_sample_size = int(sample_rate * chunk_size)
    end_index = chunk_sample_size

    audio_output = list()
    audio_fs = None

    while end_index < signal_length:
        # Bring the end index back in case it has gone over.
        end_index = min(end_index, signal_length - 1)

        signal_crop = signal[start_index:end_index]
        audio_crop, audio_fs = demodulate_signal(signal_crop, sample_rate, band_filter_cutoff, downsample_rate)
        audio_output.extend(audio_crop.tolist())

        start_index = end_index
        end_index += chunk_sample_size

    return np.array(audio_output), audio_fs
