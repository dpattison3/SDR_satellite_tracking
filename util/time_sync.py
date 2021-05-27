import numpy as np
from scipy.signal import correlate
from scipy.signal import correlation_lags


def compute_time_lag_in_samples(signal0, signal1):
    """
    Compute the time delta between two input signals by maximizing the cross correlation.
    NOTE(dominic): This is reasonably sensitive to the input signals themselves. It's often worth
    plotting the resulting correlation values.

    Args:
        signal0: first signal to sync to
        signal1: other signal to compare

    Returns:
        The lag of signal0 to signal1 in samples. When positive, signal0 is delayed from signal1.
    """
    cross_correlation = correlate(signal0, signal1, mode="full")
    lags = correlation_lags(signal0.size, signal1.size, mode="full")
    cross_correlation_lag = lags[np.argmax(cross_correlation)]

    return cross_correlation_lag


def compute_time_lag_in_seconds(signal0, signal1, sample_rate):
    """
    Compute the time delta between two input signals measured in seconds.
    See compute_time_lag_in_samples above.

    Args:
        signal0: first signal to sync to
        signal1: other signal to compare

    Returns:
        The lag of signal0 to signal1 in seconds.
    """
    time_lag_samples = compute_time_lag_in_samples(signal0, signal1)
    return time_lag_samples / sample_rate


def shift_signals_with_lag(signal0, signal1, signal0_lag):
    """
    Compute the shifted signals given the signal lag. If the lag is positive, then we zero pad
    signal1 so they align. If the lag is negative, then we zero pad signal0.

    Args:
        signal0: first signal to sync to
        signal1: other signal to compare
        signal0_lag: lag of signal0 relative to signal1 measured in samples

    Returns:
        Tuple of shifted signal0 and signal1.
    """
    lag = abs(int(signal0_lag))

    if signal0_lag >= 0:
        signal0_shifted = np.concatenate([signal0, np.zeros(lag)])
        signal1_shifted = np.concatenate([np.zeros(lag), signal1])
    else:
        signal0_shifted = np.concatenate([np.zeros(lag), signal0])
        signal1_shifted = np.concatenate([signal1, np.zeros(lag)])

    return signal0_shifted, signal1_shifted

