import numpy as np


def wrap_phase(phase):
    """
    Take a phase angle and returns a phase angle in the range [-2pi,2pi].
    """
    while phase >= (2*np.pi):
        phase -= 2*np.pi
    while phase < (-2*np.pi):
        phase += 2*np.pi
    return phase


def get_pll_filter_coefficients(frequency_bandwidth, fs):
    """
    Computes alpha and beta filter coefficients for use in a PLL.

    Args:
        frequency_bandwidth: the bandwidth within which the signal must be
        fs: sample frequency in hz

    Returns:
        tuple of alpha and beta
    """
    # A sensible default.
    damping = np.sqrt(2) / 2

    # Frequencies are in rad/sample
    loop_bw = 2*np.pi*frequency_bandwidth / fs

    coefficient_denominator = (1.0 + 2.0*damping*loop_bw + loop_bw*loop_bw)

    alpha = (4 * damping * loop_bw) / coefficient_denominator
    beta = (4 * loop_bw * loop_bw) / coefficient_denominator

    return alpha, beta


def phase_lock_loop(input_signal, fs, initial_frequency_estimate, frequency_bandwidth,
                    output_frequency_multiplier=1):
    """
    Implemented a phase lock loop.

    Heavily inspired by GNU-radio's and Lua-radio's PLL implementations.

    Args:
        input_signal: complex time series
        fs: sample rate in hz
        initial_frequency_estimate: initial guess of the frequency to lock to
        frequency_bandwidth: how far the PLL is expected to stray from the initial guess
        output_frequency_multiplier: creates a harmonic that is in-phase to the input

    Returns:
        tuple of:
            the output signal that is frequency and phase matched to the input
            a harmonic of that output dicated by output_frequency_multiplier
            the error signal (could be used for FM demodulation)
    """
    phase = 0.0
    phase_harmonic = 0.0
    frequency = 2*np.pi * initial_frequency_estimate / fs

    alpha,beta = get_pll_filter_coefficients(frequency_bandwidth, fs)

    output = np.zeros(len(input_signal), dtype=np.complex64)
    output_harmonic = np.zeros(len(input_signal), dtype=np.complex64)
    error = np.zeros(len(input_signal), dtype=np.float64)

    for i in range(len(input_signal)):
        # Output signals are the references we are trying to track to the correct frequency and
        # phase. They correspond to the voltage controlled oscillator in EE literature.
        output[i] = np.cos(phase) + 1j * np.sin(phase)
        output_harmonic[i] = np.cos(phase_harmonic) + 1j * np.sin(phase_harmonic)

        # Phase detector is differencing the phase of the input and output.
        error[i] = np.angle(input_signal[i] * np.conj(output[i]))

        # Loop filtering updates the phase and frequency of the reference. The phase always gets
        # the frequency added because the frequency is in rad/sample.
        frequency += beta * error[i]
        phase += frequency + alpha * error[i]

        # Same process for the harmonic, but we take a multiple of the increment.
        # NOTE(dominic): Lua radio seems to think that it should only be the frequency increment
        # that gets multiplied by the harmonic. That makes more sense in theory, but in practice,
        # this leaves a small constant phase offset.
        phase_harmonic += output_frequency_multiplier * (frequency + alpha * error[i])

        # wrap phi angles
        phase = wrap_phase(phase)
        phase_harmonic = wrap_phase(phase_harmonic)

    return output, output_harmonic, error
