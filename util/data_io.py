import numpy as np


def read_raw_iq_data(input_path, dtype_string, offset_value=0, start_index_fraction=None, end_index_fraction=None):
    """
    Read the raw data file output by the RTL-SDR command line tool or as recorded by GQRX.
    File is expected to be IQ data channels interleaved formated as dicated by dtype_string.
    The rtl-sdr command line tool reads int8 while GQRX reads float32.
    raw_data = [d0, d1, d2, d3, ...]
    complex_data = [(d0, i * d1), (d2, i * d3), ...]

    Optionally pass in indices to start and stop at for the original data. This allows keeping some
    fraction of the raw data to save on memory. Done using fraction of original data length so the
    user does not need to know the length before calling.
    Note, because the raw data is uint8 and the output is floats, keeping only a fraction of the
    read data can significantly increase the amount read. However, we are still reading the full
    raw data. Fixing that would require writing our own function to read the raw byte stream.

    Example usages:
    For recorded by rtl-sdr:
    complex_signal = read_raw_iq_data(input_path, dtype_string="uint8", offset_value=-127.5)

    For recorded by GQRX:
    complex_signal = read_raw_iq_data(input_path, dtype_string="float32")

    Args:
        input_path: full path to the data file
        dtype_string: the type of the encoded data, ex: uint8 or float32
        offset_value: an optional offset to apply to all data - useful for offsetting unisgned
                      values back to zero-centered
        start_index_fraction: return data will start from index floor(raw_length*start_index_fraction)
        end_index_fraction: return data will stop before index floor(raw_length*end_index_fraction)

    Returns:
        np.ndarry of complex signal as float32 values
    """
    raw_signal = np.fromfile(input_path, dtype=dtype_string)

    if start_index_fraction and (start_index_fraction < 0 or start_index_fraction >= 1):
        raise ValueError("start_index_fraction must be positive and less than 1.")
    elif end_index_fraction and (end_index_fraction < 0 or end_index_fraction >= 1):
        raise ValueError("end_index_fraction must be positive and less than 1.")
    elif start_index_fraction and end_index_fraction and \
            (end_index_fraction <= start_index_fraction):
        raise ValueError("end_index_fraction must be greater than start_index_fraction.")

    if len(raw_signal) % 2 == 1:
        raise ValueError("Invalid data file. File must be IQ with even length to represent real and"
                         " imaginary data pairs.")

    raw_signal_length = len(raw_signal) / 2
    start_index = None
    end_index = None

    print(f"loaded raw signal of length {int(raw_signal_length)}")

    if start_index_fraction:
        start_index = 2 * np.floor(raw_signal_length * start_index_fraction).astype(np.int)

    if end_index_fraction:
        end_index = 2 * np.floor(raw_signal_length * end_index_fraction).astype(np.int)

    merged_signal = raw_signal[start_index:end_index].astype(np.float32) - offset_value
    real_signal = merged_signal[::2]
    imag_signal = merged_signal[1::2]

    return real_signal + 1j * imag_signal


def read_rtl_raw_data(input_path, start_index_fraction=None, end_index_fraction=None):
    """
    Read the raw data file output by the RTL-SDR tool. This is assumed to be encoded as unsigned
    8-bit integers. We offset by -127.5 to get the signal zero-centered.
    """
    return read_raw_iq_data(
        input_path, dtype_string="uint8", offset_value=-127.5,
        start_index_fraction=start_index_fraction, end_index_fraction=end_index_fraction
    )


def read_gqrx_raw_data(input_path, start_index_fraction=None, end_index_fraction=None):
    """
    Read the raw data file output by the recording GQRX. This is assumed to be encoded as 32-bit
    floating point values. Since these are signed values, we do not need to offset them.
    """
    return read_raw_iq_data(
        input_path, dtype_string="float32", start_index_fraction=start_index_fraction,
        end_index_fraction=end_index_fraction
    )
