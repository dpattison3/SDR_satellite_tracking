from scipy.fftpack import fft, fftfreq, next_fast_len
from scipy.signal import sosfreqz
import numpy as np

from util.filtering import compute_lpf_coeff


def compute_fft_plot(signal, sampling_dt, center_frequency=0.0, pad=True):
    """
    Compute the FFT of a complex signal and its corresponding frequencies used.
    Assumes signal is complex and is a 1-D array.
    Optionally, shift the frequency.

    Args:
        signal: numpy.ndarray of 1-D complex signal to use
        sampling_dt: float of inverse of the sample rate
        center_frequency: float representing the center frequency of a demodulated signal
        pad: when true, zero pad to the next 5-smooth number
             see scipy's documentation on next_fast_len

    Returns:
        tuple of numpy.ndarray with the first being
    """

    # Optionally, pad the size. This can greatly increase the speed of the FFT algorithm, several
    # orders of magnitude in the worst cases.
    N = len(signal)
    if pad:
        N = next_fast_len(N)

    # Compute the FFT magnitudes and scale them appropriately.
    fft_frequency_magnitudes = fft(signal, n=N)
    fft_frequency_magnitudes = 2.0 / N * np.abs(fft_frequency_magnitudes)
    # Compute the corresponding frequencies.
    fft_frequencies = fftfreq(N, sampling_dt) + center_frequency

    # Shift the frequencies and their corresponding magnitudes because, by default, they start from
    # 0, list positive values, then list negative values. The shift makes the frequencies
    # monotonically increasing, perfect for plotting.
    return np.fft.fftshift(fft_frequencies), np.fft.fftshift(fft_frequency_magnitudes)


def compute_fft_plot_from_sample_rate(signal, sampling_rate, center_frequency=0.0, pad=True):
    """
    see compute_fft_plot above

    Args:
        signal: numpy.ndarray of 1-D complex signal to use
        sampling_rate: float or int of the sample rate
        center_frequency: float representing the center frequency of a demodulated signal
        pad: when true, allow zero padding to speed up FFT computation

    Returns:
        tuple of numpy.ndarray with the first being
    """
    return compute_fft_plot(signal, 1.0 / sampling_rate, center_frequency=center_frequency, pad=pad)


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

