#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <cufft.h>         // Include cuFFT header for cufftComplex
#include <cuda_runtime.h>  // Include CUDA runtime header
#include "Array.h"

extern "C" void fft3d(float* h_data, int nx, int ny, int nz, bool inverse);

void print_values(const float* data, int nx, int ny, int nz) {
    for (int i = 0; i < 10; ++i) {
        std::cout << data[i] << " ";
    }
    std::cout << std::endl;
}

void compare_values(const float* original, const float* transformed, int nx, int ny, int nz) {
    for (int i = 0; i < 10; ++i) {
        std::cout << "Original: " << original[i] << " Transformed: " << transformed[i] << std::endl;
    }
}

int main() {
    const int nx = 256;
    const int ny = 256;
    const int nz = 256;
    size_t size = nx * ny * nz * sizeof(float);
    size_t complex_size = nx * ny * (nz / 2 + 1) * sizeof(cufftComplex);
    size_t nzp=nz/2+1;

    Array::Array3<cufftComplex> data_c(nx,ny,nzp);
    Array::Array3<cufftReal>    original_data_c(nx,ny,nz);
    cufftReal * data = reinterpret_cast<cufftReal *>(&data_c[0][0][0]);
    cufftReal * original_data = &original_data_c[0][0][0];

    // float* data = (float*)malloc(complex_size);
    //float* original_data = (float*)malloc(size);
    fft3d(data, nx, ny, nz, false);
    

    // Initialize the data
    for (int i = 0; i < nx * ny * nz; ++i) {
        data[i] = static_cast<float>(rand()) / RAND_MAX;
        original_data[i] = data[i];
    }

    // Timing the forward FFT
    auto start = std::chrono::high_resolution_clock::now();
    fft3d(data, nx, ny, nz, false);
    

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double forward_time = duration.count();
    std::cout << "Forward FFT time: " << forward_time << " seconds" << std::endl;

    // Timing the inverse FFT
    start = std::chrono::high_resolution_clock::now();
    fft3d(data, nx, ny, nz, true);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    double inverse_time = duration.count();
    std::cout << "Inverse FFT time: " << inverse_time << " seconds" << std::endl;

    // Calculate performance in MFlop/second
    // An FFT has O(N log N) operations
    // Here, N = nx * ny * nz
    size_t N = nx * ny * nz;
    double num_ops = 5.0 * N * log2(N); // FFT flop count approximation
    double forward_performance = num_ops / (forward_time * 1e6);
    double inverse_performance = num_ops / (inverse_time * 1e6);

    std::cout << "Forward FFT performance: " << forward_performance << " MFlop/s" << std::endl;
    std::cout << "Inverse FFT performance: " << inverse_performance << " MFlop/s" << std::endl;

    // Compare values
    compare_values(original_data, data, nx, ny, nz);

    // Clean up
    // free(data);
    // free(original_data);

    return 0;
}
