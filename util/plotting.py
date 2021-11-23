from scipy.fftpack import fft, fftfreq
from scipy.signal import sosfreqz
import numpy as np

from util.filtering import compute_lpf_coeff


def compute_fft_plot(signal, sampling_dt, center_frequency=0.0):
    """
    Compute the FFT of a complex signal and its corresponding frequencies used.
    Assumes signal is complex and is a 1-D array.
    Optionally, shift the frequency.

    Args:
        signal: numpy.ndarray of 1-D complex signal to use
        sampling_dt: float of inverse of the sample rate
        center_frequency: float representing the center frequency of a demodulated signal

    Returns:
        tuple of numpy.ndarray with the first being
    """

    num_samples = len(signal)
    # Compute the FFT magnitudes and scale them appropriately.
    fft_frequency_magnitudes = fft(signal)
    fft_frequency_magnitudes = 2.0 / num_samples * np.abs(fft_frequency_magnitudes)
    # Compute the corresponding frequencies.
    fft_frequencies = fftfreq(num_samples, sampling_dt) + center_frequency

    # Shift the frequencies and their corresponding magnitudes because, by default, they start from
    # 0, list positive values, then list negative values. The shift makes the frequencies
    # monotonically increasing, perfect for plotting.
    return np.fft.fftshift(fft_frequencies), np.fft.fftshift(fft_frequency_magnitudes)


def compute_fft_plot_from_sample_rate(signal, sampling_rate, center_frequency=0.0):
    """
    see compute_fft_plot above

    Args:
        signal: numpy.ndarray of 1-D complex signal to use
        sampling_rate: float or int of the sample rate
        center_frequency: float representing the center frequency of a demodulated signal

    Returns:
        tuple of numpy.ndarray with the first being
    """
    return compute_fft_plot(signal, 1.0 / sampling_rate, center_frequency=center_frequency)


def compute_frequency_response(sos, sample_rate):
    """
    Computes the frequency response of a given filter with sos coefficients.

    Args:
        sos: filter parameters
        sample_rate: sample rate of the given filter

    Returns:
        tuple of frequencies and their respective magnitudes in response to the filter
    """
    w, h = sosfreqz(sos, worN=2000)

    frequencies = (sample_rate * 0.5 / np.pi) * w
    magnitudes = abs(h)

    return frequencies, magnitudes


def compute_lpf_frequency_response(cutoff_frequency, sample_rate, order=10):
    """
    Computes the frequency response of a low pass filter with parameters as given.


    Args:
        cutoff_frequency: frequency of the cutoff of the LPF
        sample_rate: sample rate of the given filter
        order: order of the filter

    Returns:
        tuple of frequencies and their respective magnitudes in response to the filter
    """

    sos = compute_lpf_coeff(cutoff_frequency, sample_rate, order)
    return compute_frequency_response(sos, sample_rate)

