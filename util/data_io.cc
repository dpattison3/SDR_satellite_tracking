#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>


extern "C" {
    int ReadRawFile(const char* const file_path, const int size, const float offset_frequency,
        float sampling_rate, const float low_pass_cutoff_frequency, const int downsample_factor,
        float* const real_output, float* const imaginary_output);
}


class LowPassFilter {
 public:
  static float ComputeHammingCoefficient(const int n, const int length) {
    return 0.54 - 0.46 * std::cos((2*M_PI*n) / (length-1));
  }

  // The number of coefficients requested must be even!
  LowPassFilter(const int num_coefficients, const float cutoff_frequency, const float sampling_frequency)
      : num_coefficients_(num_coefficients) {
    // Cutoff is normalized by the nyquist rate.
    const float normalized_cutoff_frequency = cutoff_frequency / (sampling_frequency / 2);

    coefficients_.reserve(num_coefficients);

    float raw_coefficient_sum = 0;
    const float ntaps = static_cast<float>(num_coefficients);
    for (int i = 0; i < num_coefficients; ++i) {
      const float coeff = std::sin( M_PI * normalized_cutoff_frequency * (i - (ntaps-1)/2) ) / (M_PI * (i - (ntaps-1)/2));
      const float hamming_coeff = ComputeHammingCoefficient(i, num_coefficients);

      const float coefficient_scaled = coeff * hamming_coeff;
      coefficients_.push_back(coefficient_scaled);
      raw_coefficient_sum += coefficient_scaled;
    }

    for (float& coeff : coefficients_) {
        coeff /= raw_coefficient_sum;
    }

    input_buffer_.resize(num_coefficients);
    std::fill(input_buffer_.begin(), input_buffer_.end(), std::complex<float>{0, 0});
  }

  std::complex<float> Update(const std::complex<float> input_value) {
    input_buffer_[0] = input_value;

    std::complex<float> output_value = {0, 0};
    for (int i = 0; i < num_coefficients_; ++i) {
      output_value += input_buffer_[i] * coefficients_[i];
    }

    // Now we can shift the buffer.
    // NOTE(dominic): This is super wasteful. We could do better if we used a ring buffer, but
    // I can't be bothered...
    // We go backwards to avoid propagating one value downwards.
    for (int i = (num_coefficients_-1); i > 0; --i) {
        input_buffer_[i] = input_buffer_[i-1];
    }

    return output_value;
  }

 private:
  const int num_coefficients_;
  std::vector<float> coefficients_;
  std::vector<std::complex<float>> input_buffer_;

};

/*
def compute_hamming_coefficient(n, M):
    return 0.54 - 0.46*np.cos((2*np.pi*n)/(M-1))


def compute_fir_coefficients(num_taps, cutoff_frequency, sampling_frequency):
    if (num_taps % 2) != 0:
        raise RuntimeError("Number of coefficients in this FIR must be even")

    nyquist_frequency = sampling_frequency / 2
    normalized_cutoff = cutoff_frequency / nyquist_frequency

    coefficients = []
    for i in range(num_taps):
        coefficient = np.sin(np.pi*normalized_cutoff*(i - (num_taps-1)/2))/(np.pi*(i - (num_taps-1)/2))
        hamming_coefficient = compute_hamming_coefficient(i, num_taps)

        coefficients.append(-1 * coefficient * hamming_coefficient)

    total_sum = np.sum(coefficients)
    coefficients = [c / total_sum for c in coefficients]

    return coefficients
*/


int ReadRawFile(
        const char* const file_path, const int size, const float offset_frequency,
        float sampling_rate, const float low_pass_cutoff_frequency, const int downsample_factor,
        float* const real_output, float* const imaginary_output) {

    // A handy constant.
    constexpr std::complex<float> double_i_pi = {0, 2*M_PI};

    printf("Reading from file path %s.\nWe will read %i pairs of floats.\nLowpass filter cutoff is %f\n", file_path, size, low_pass_cutoff_frequency);
    std::ifstream fin(file_path, std::ios::binary);

    LowPassFilter lpf{/* num_coefficients */ 100, low_pass_cutoff_frequency, sampling_rate};

    float real_value, imaginary_value;
    std::complex<float> input_value, shift_value, shifted_input_value;
    int num_values_read = 0, num_values_written = 0;

    for (int i = 0; i < size; ++i) {
        // Read a pair of values representing the real and imaginary component. If either fail
        // because we have reached the end of the file before the desired number of values to read,
        // then the user has goofed. Warn them!
        if (!fin.read(reinterpret_cast<char*>(&real_value), sizeof(float))) {
            printf(
                "Requested to read %i pairs of floats, but we reached the end of the file first.\n",
                size);
            break;
        }
        if (!fin.read(reinterpret_cast<char*>(&imaginary_value), sizeof(float))) {
            printf(
                "Requested to read %i pairs of floats, but we reached the end of the file first.\n",
                size);
            break;
        }
        ++num_values_read;

        input_value = {real_value, imaginary_value};
        shift_value = std::exp(double_i_pi * (-offset_frequency * i / sampling_rate));
        shifted_input_value = input_value * shift_value;

        const std::complex<float> lpf_output = lpf.Update(shifted_input_value);
        // printf("output %f  %f\n", lpf_output.real(), lpf_output.imag());

        if (i % downsample_factor == 0) {
          real_output[num_values_written] = lpf_output.real();
          imaginary_output[num_values_written] = lpf_output.imag();
          ++num_values_written;
        }
    }


    printf("Finished reading file. Read %d and wrote %d.\n", num_values_read, num_values_written);

    return num_values_read;
}
