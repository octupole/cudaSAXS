#include <iostream>
#include <cufft.h>
#include <cuda_runtime.h>

extern "C" {
#define CHECK_CUDA_ERROR(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err) << " at line " << __LINE__ << std::endl; \
            exit(err); \
        } \
    } while (0)

#define CHECK_CUFFT_ERROR(call) \
    do { \
        cufftResult err = call; \
        if (err != CUFFT_SUCCESS) { \
            std::cerr << "CUFFT Error: " << err << " at line " << __LINE__ << std::endl; \
            exit(err); \
        } \
    } while (0)

void fft3d(float* h_data, int nx, int ny, int nz, bool inverse) {
    cufftHandle plan;
    cufftComplex* d_data;
    size_t size = nx * ny * nz * sizeof(float);
    size_t complex_size = nx * ny * (nz / 2 + 1) * sizeof(cufftComplex);

    // Allocate device memory
    CHECK_CUDA_ERROR(cudaMalloc((void**)&d_data, complex_size));
    if (inverse) {
        // Copy data to device
        CHECK_CUDA_ERROR(cudaMemcpy(d_data, h_data, complex_size, cudaMemcpyHostToDevice));

        // Create a 3-D FFT plan for inverse transform
        CHECK_CUFFT_ERROR(cufftPlan3d(&plan, nx, ny, nz, CUFFT_C2R));

        // Execute the inverse FFT
        CHECK_CUFFT_ERROR(cufftExecC2R(plan, d_data, (cufftReal*)d_data));

        // Copy the result back to host
        CHECK_CUDA_ERROR(cudaMemcpy(h_data, d_data, size, cudaMemcpyDeviceToHost));

        // Normalize the result
        for (int i = 0; i < nx * ny * nz; ++i) {
            h_data[i] /= (nx * ny * nz);
        }
    } else {
        // Copy data to device
        CHECK_CUDA_ERROR(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice));
        // Create a 3-D FFT plan for forward transform
        CHECK_CUFFT_ERROR(cufftPlan3d(&plan, nx, ny, nz, CUFFT_R2C));

        // Execute the forward FFT
        CHECK_CUFFT_ERROR(cufftExecR2C(plan, (cufftReal*)d_data, d_data));


        // Copy the result back to host
        CHECK_CUDA_ERROR(cudaMemcpy(h_data, d_data, complex_size, cudaMemcpyDeviceToHost));
        }

    // Clean up
    cufftDestroy(plan);
    cudaFree(d_data);
}
}