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
    sos = butter(10, cutoff_frequency / nyquist_rate, btype='low', output='sos')
    return sos


def filter_complex_signal(signal, cutoff_frequency, sample_rate, order=10):
    """
    Compute the low pass filtered signal of the complex input signal.
    This splits the signal into a real and imaginary component, filters them separately, then
    rejoins them into a complex filtered output.

    Args:
        signal: complex signal to filter
        cutoff_frequency: the frequency of the LPF cutoff
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter

    Returns:
        complex output signal that has been low pass filtered
    """

    # We have a complex signal but a real filter. Split the data, filter separately, then re-combine.
    
    mean_shifted_signal = signal - np.mean(signal)
    real_signal = np.real(mean_shifted_signal)
    imag_signal = np.imag(mean_shifted_signal)
    
    sos = compute_lpf_coeff(cutoff_frequency, sample_rate, order=order)

    filtered_real = sosfilt(sos, real_signal)
    filtered_imag = sosfilt(sos, imag_signal)

    filtered_signal = filtered_real + 1j * filtered_imag
    return filtered_signal


def filter_real_signal(signal, cutoff_frequency, sample_rate, order=10):
    """
    Compute the low pass filtered signal of the real input signal.
    see filter_complex_signal above.

    Args:
        signal: complex signal to filter
        cutoff_frequency: the frequency of the LPF cutoff
        sample_rate: expected sample rate of the signal to be filtered
        order: order of the filter

    Returns:
        real output signal that has been low pass filtered
    """

    # NOTE(dominic): Technically, this is a bit more expensive, but it's clean and I suspect for now
    # it's fine.
    filtered_signal = filter_complex_signal(signal, cutoff_frequency, sample_rate, order)
    return np.real(filtered_signal)
