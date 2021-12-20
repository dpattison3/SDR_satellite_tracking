import numpy as np
from scipy.signal import butter, sosfilt


def compute_lpf_coeff(cutoff_frequency, sample_rate, order=10):
    """
    Compute the LPF sos coefficients of an LPF using scipy's butter function.

    Args:
        cutoff_frequency: the frequency of the LPF cutoff
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter

    Returns:
        sos coefficients of the filter
    """
    nyquist_rate = 0.5 * sample_rate
    sos = butter(order, cutoff_frequency / nyquist_rate, btype='low', output='sos')
    return sos


def compute_bandpass_coeff(low_cutoff, high_cutoff, sample_rate, order=10):
    """
    Compute coefficients for a bandpass filter using scipy's butter function.

    Args:
        low_cutoff: the frequency that starts the band pass
        high_cutoff: the frequency that ends the band pass
        order: order of the filter

    Returns:
        sos coefficients of the filter
    """
    nyquist_rate = 0.5 * sample_rate
    start_frequency = low_cutoff / nyquist_rate
    stop_frequency = high_cutoff / nyquist_rate
    sos = butter(order, [start_frequency, stop_frequency], btype='bandpass', output='sos')
    return sos


def filter_complex_signal(signal, filter_sos_coefficients, apply_mean_shift=True):
    """
    Compute the filtered signal of the complex input signal.
    This splits the signal into a real and imaginary component, filters them separately, then
    rejoins them into a complex filtered output.

    Args:
        signal: complex signal to filter
        filter_sos_coefficients: filter coefficients to pass to scipy's sosfilt
        apply_mean_shift: remove the DC component of the raw signal before filtering

    Returns:
        complex output signal that has been filtered
    """
    mean_shifted_signal = (signal - np.mean(signal)) if apply_mean_shift else signal

    filtered_real = sosfilt(filter_sos_coefficients, np.real(mean_shifted_signal))
    filtered_imag = sosfilt(filter_sos_coefficients, np.imag(mean_shifted_signal))

    return filtered_real + 1j * filtered_imag


def low_pass_filter_complex_signal(signal, cutoff_frequency, sample_rate, order=10,
                                   apply_mean_shift=True):
    """
    Compute the low pass filtered signal of the complex input signal.
    This splits the signal into a real and imaginary component, filters them separately, then
    rejoins them into a complex filtered output.

    Args:
        signal: complex signal to filter
        cutoff_frequency: the frequency of the LPF cutoff
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter
        apply_mean_shift: remove the DC component of the raw signal before filtering

    Returns:
        complex output signal that has been low pass filtered
    """

    sos = compute_bandpass_coeff(cutoff_frequency, sample_rate, order=order)
    return filter_complex_signal(signal, sos, apply_mean_shift)


def band_pass_filter_complex_signal(signal, low_cutoff_frequency, high_cutoff_frequency,
                                    sample_rate, order=10, apply_mean_shift=True):
    """
    Compute the low pass filtered signal of the complex input signal.
    This splits the signal into a real and imaginary component, filters them separately, then
    rejoins them into a complex filtered output.

    Args:
        signal: complex signal to filter
        low_cutoff_frequency: the frequency of the start of the bandpass
        high_cutoff_frequency: the frequency of the start of the bandpass
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter
        apply_mean_shift: remove the DC component of the raw signal before filtering

    Returns:
        complex output signal that has been low pass filtered
    """
    sos = compute_bandpass_coeff(low_cutoff_frequency, high_cutoff_frequency, sample_rate,
                                 order=order)
    return filter_complex_signal(signal, sos, apply_mean_shift)


def low_pass_filter_real_signal(signal, cutoff_frequency, sample_rate, order=10,
                                apply_mean_shift=True):
    """
    Compute the low pass filtered signal of the real input signal.
    see filter_complex_signal above.

    Args:
        signal: complex signal to filter
        cutoff_frequency: the frequency of the LPF cutoff
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter
        apply_mean_shift: remove the DC component of the raw signal before filtering

    Returns:
        real output signal that has been low pass filtered
    """

    # NOTE(dominic): Technically, this is a bit more expensive, but it's clean and I suspect for now
    # it's fine.
    filtered_signal = low_pass_filter_complex_signal(signal, cutoff_frequency, sample_rate, order,
                                                     apply_mean_shift)
    return np.real(filtered_signal)
