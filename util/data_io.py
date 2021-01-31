import numpy as np


def read_sdriq_data(input_path):
    """
    Read the raw data file output by the RTL-SDR tool.
    File is expected to be IQ data formated as int8 alternating between real and associated
    imaginary components.
    raw_data = [d0, d1, d2, d3, ...]
    complex_data = [(d0, i * d1), (d2, i * d3), ...]

    Args:
        input_path: full path to the data file
    
    Returns:
        np.ndarry of complex signal
    """
    raw_signal = np.fromfile(input_path, dtype="uint8")

    merged_signal = raw_signal - 127.5
    real_signal = merged_signal[::2]
    imag_signal = merged_signal[1::2]

    return real_signal + 1j * imag_signal
