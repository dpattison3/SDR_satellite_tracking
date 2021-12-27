#include <stdio.h>
#include <math.h>
#include <complex.h>


#define PI2 ((float)(M_PI * 2))


void ComputeControlLoopCoefficients(const float frequency_bandwidth, const float sampling_rate,
                                    float* const Kp, float* const Ki) {
  // TODO(dominic): Expose this as a parameter. For now, I think this serves our usecase. It
  // represents the damping coefficient of the system.
  const float zeta = 1/ sqrt(2);

  // K0 can just be fixed to 1. It represents an arbitrary gain factor.
  // Kd depends on the gain that needs to be applied for a given phase error computation. In our
  // case, this is just 1.
  const float Kd = 1.0f;
  const float K0 = 1.0f;

  const float normalized_frequency_bandwidth = frequency_bandwidth / sampling_rate;

  const float denominator = (zeta + (1/(4*zeta)));

  *Kp = 1/(Kd*K0) * 4*zeta / denominator * normalized_frequency_bandwidth;
  *Ki = 1/(Kd*K0) * 4 / (denominator*denominator) * \
      normalized_frequency_bandwidth * normalized_frequency_bandwidth;
}

int PhaseLockedLoop(size_t num_samples, const float* const real_input, const float* const imaginary_input,
                    float* const real_output, float* const imaginary_output,
                    float* const output_frequencies, const float sampling_rate,
                    const float initial_frequency_estimate, const float frequency_bandwidth) {
  float Kp=0,Ki=0;
  ComputeControlLoopCoefficients(frequency_bandwidth, sampling_rate, &Kp, &Ki);

  float frequency_estimate = PI2 * initial_frequency_estimate / sampling_rate;
  float phase_estimate = 0.0f, phase_error = 0.0f;
  float complex complex_input, complex_output;
  size_t num_samples_processed = 0;

  for (int i = 0; i < num_samples; ++i) {
    // Use the current phase estimate to update the output signal.
    real_output[i] = cos(phase_estimate);
    imaginary_output[i] = sin(phase_estimate);

    // Compute the phase error between the output and input. This is just the phase angle
    // difference. We use imaginary numbers here for convenience.
    complex_input = real_input[i] + I * imaginary_input[i];
    complex_output = real_output[i] + I * imaginary_output[i];
    phase_error = cargf(complex_input * conjf(complex_output));

    // Loop filtering updates the phase and frequency of the reference. The phase always gets
    // the frequency added because the frequency is in rad/sample.
    frequency_estimate += Ki * phase_error;
    output_frequencies[i] = frequency_estimate;

    phase_estimate += frequency_estimate + Kp * phase_error;

    ++num_samples_processed;
  }

  return num_samples_processed;
}

// NOTE(dominic): This seems a little repetitive, but we want this to be as fast as possible, so
// skipping the harmonic computation altogether is better than having a conditional or wasted
// computation in the processing loop.
int PhaseLockedLoopWithHarmonic(
      size_t num_samples, const float* const real_input, const float* const imaginary_input,
      float* const real_output, float* const imaginary_output, float* const real_output_harmonic,
      float* const imaginary_output_harmonic, float* const output_frequencies,
      const float sampling_rate, const float initial_frequency_estimate,
      const float frequency_bandwidth, const int harmonic_multiple) {
  float Kp=0,Ki=0;
  ComputeControlLoopCoefficients(frequency_bandwidth, sampling_rate, &Kp, &Ki);

  float frequency_estimate = PI2 * initial_frequency_estimate / sampling_rate;
  float phase_estimate = 0.0f, phase_harmonic_estimate = 0.0f, phase_error = 0.0f,
        phase_increment = 0.0f;
  float complex complex_input, complex_output;
  size_t num_samples_processed = 0;

  for (int i = 0; i < num_samples; ++i) {
    // Use the current phase estimate to update the output signal.
    real_output[i] = cos(phase_estimate);
    imaginary_output[i] = sin(phase_estimate);

    real_output_harmonic[i] = cos(phase_harmonic_estimate);
    imaginary_output_harmonic[i] = sin(phase_harmonic_estimate);

    // Compute the phase error between the output and input. This is just the phase angle
    // difference. We use imaginary numbers here for convenience.
    complex_input = real_input[i] + I * imaginary_input[i];
    complex_output = real_output[i] + I * imaginary_output[i];
    phase_error = cargf(complex_input * conjf(complex_output));

    // Loop filtering updates the phase and frequency of the reference. The phase always gets
    // the frequency added because the frequency is in rad/sample.
    frequency_estimate += Ki * phase_error;
    output_frequencies[i] = frequency_estimate;

    phase_increment = frequency_estimate + Kp * phase_error;
    phase_estimate += phase_increment;
    // NOTE(dominic): This doesn't make the most intuitive sense as we only actually want to double
    // the frequency increment. Empirically, this produces a result that is not offset in phase.
    phase_harmonic_estimate += harmonic_multiple * phase_increment;

    ++num_samples_processed;
  }

  return num_samples_processed;
}
