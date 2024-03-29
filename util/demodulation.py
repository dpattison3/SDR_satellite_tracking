import numpy as np

from util.filtering import low_pass_filter_complex_signal, low_pass_filter_real_signal


# Cutoff used in filtering the audio signal. Since this is the top of the audio spectrum, we just
# define it as a constant.
AUDIO_LPF_FREQUENCY = 20e3


def demodulate_signal(signal, sample_rate, base_band_filter_cutoff=60000,
                      base_band_downsample_rate=10, audio_filter_cutoff=15000,
                      audio_downsample_rate=5,
                      apply_output_filter=True):
    """
    Demodulates an FM signal. First, perform a LPF of input to the appropriate bandwidth. Next,
    diff the phase angles to get the data signal. LPF that signal to just the audible range.
    Finally, downsample to a more reasonable sample rate.

    Args:
        signal: complex baseband input signal
        sample_rate: sample rate of raw signal
        base_band_filter_cutoff: the first LPF of raw IQ signal
        base_band_downsample_rate: downsampling of the raw IQ signal after filtering
        audio_filter_cutoff: LPF of the demodulated signal
        audio_downsample_rate: downsample the demodulated signal after filtering
        apply_output_filter: when true, apply a LPF at audio_filter_cutoff to the demodulated signal

    Returns:
        tuple of output audio signal and its associated sample rate
    """

    filtered_signal = low_pass_filter_complex_signal(
        signal=signal, cutoff_frequency=base_band_filter_cutoff, sample_rate=sample_rate)

    filtered_signal_downsampled = filtered_signal[::base_band_downsample_rate]
    filtered_signal_sample_rate = int(sample_rate / base_band_downsample_rate)

    # Compute the phase angle difference between temporally successive values in the complex signal.
    # We use the product of the complex conjugate because that efficiently handles the wraparound
    # that you would get if you just differenced phase angles.
    angle_diff = np.angle(
        np.conjugate(filtered_signal_downsampled[:-1]) * filtered_signal_downsampled[1:]
    )

    filtered_phase_angle_diff = angle_diff
    if apply_output_filter:
        filtered_phase_angle_diff = low_pass_filter_real_signal(
            signal=angle_diff, cutoff_frequency=audio_filter_cutoff,
            sample_rate=filtered_signal_sample_rate
        )

    phase_angle_diff_downsampled = filtered_phase_angle_diff[::audio_downsample_rate]
    phase_angle_diff_downsampled_sample_rate = \
        int(filtered_signal_sample_rate / audio_downsample_rate)

    return phase_angle_diff_downsampled, phase_angle_diff_downsampled_sample_rate


def chunked_demodulate_signal(signal, sample_rate, chunk_size=1, base_band_filter_cutoff=60000,
                              base_band_downsample_rate=10, audio_filter_cutoff=15000,
                              audio_downsample_rate=5,
                              apply_output_filter=True):
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
        base_band_filter_cutoff: the first LPF of raw IQ signal
        base_band_downsample_rate: downsampling of the raw IQ signal after filtering
        audio_filter_cutoff: LPF of the demodulated signal
        audio_downsample_rate: downsample the demodulated signal after filtering
        apply_output_filter: when true, apply a LPF at audio_filter_cutoff to the demodulated signal

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
        audio_crop, audio_fs = demodulate_signal(signal_crop, sample_rate, base_band_filter_cutoff,
                                                 base_band_downsample_rate, audio_filter_cutoff,
                                                 audio_downsample_rate, apply_output_filter)
        audio_output.extend(audio_crop.tolist())

        start_index = end_index
        end_index += chunk_sample_size

    return np.array(audio_output), audio_fs
