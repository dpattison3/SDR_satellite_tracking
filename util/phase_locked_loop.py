import ctypes
import numpy as np


# This path must be relative to where the python code is executed from.
c_phase_locked_loop = ctypes.CDLL('build/phase_locked_loop.so')


def cast_list_to_cfloat(input: np.ndarray):
    contiguous_array = np.ascontiguousarray(input, dtype=np.float32)
    return (ctypes.c_float * len(contiguous_array)).from_buffer(contiguous_array)


def phase_locked_loop(input_signal, sampling_rate, initial_frequency_estimate,
                      frequency_bandwidth, output_frequency_harmonic=None):
    real_signal = np.real(input_signal)

    real_signal = cast_list_to_cfloat(real_signal)
    imaginary_signal = cast_list_to_cfloat(np.imag(input_signal))

    real_output = cast_list_to_cfloat(np.zeros_like(real_signal))
    imaginary_output = cast_list_to_cfloat(np.zeros_like(real_signal))

    frequencies = cast_list_to_cfloat(np.zeros_like(real_signal))

    if output_frequency_harmonic:
        real_output_harmonic = cast_list_to_cfloat(np.zeros_like(real_signal))
        imaginary_output_harmonic = cast_list_to_cfloat(np.zeros_like(real_signal))

        c_phase_locked_loop.PhaseLockedLoopWithHarmonic(
            len(real_signal), real_signal, imaginary_signal, real_output, imaginary_output,
            real_output_harmonic, imaginary_output_harmonic, frequencies,
            ctypes.c_float(sampling_rate), ctypes.c_float(initial_frequency_estimate),
            ctypes.c_float(frequency_bandwidth),output_frequency_harmonic
        )
    else:
        c_phase_locked_loop.PhaseLockedLoop(
            len(real_signal), real_signal, imaginary_signal, real_output, imaginary_output,
            frequencies, ctypes.c_float(sampling_rate), ctypes.c_float(initial_frequency_estimate),
            ctypes.c_float(frequency_bandwidth)
        )

    output = np.array(real_output) + 1j * np.array(imaginary_output)
    frequencies = np.array(frequencies)

    if output_frequency_harmonic:
        output_harmonic = np.array(real_output_harmonic) + 1j * np.array(imaginary_output_harmonic)
    else:
        output_harmonic = None

    return output, output_harmonic, frequencies