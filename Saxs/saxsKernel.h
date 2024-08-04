/**
 * @class saxsKernel
 * @brief Provides functionality for performing small-angle X-ray scattering (SAXS) calculations.
 *
 * The `saxsKernel` class is responsible for managing the memory and computation required for SAXS
 * calculations. It provides methods for setting the number of particles and grid dimensions, as well
 * as running the main SAXS kernel. The class also manages the allocation and deallocation of
 * various device memory buffers used in the SAXS computations.
 */
#ifndef SAXSKERNEL_H
#define SAXSKERNEL_H
#include "Splines.h"
#include "Options.h"
#include <vector>
#include <cufft.h>
#include <cuComplex.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>
#include <map>
#pragma once

class saxsKernel
{
public:
    saxsKernel(int _nx, int _ny, int _nz, int _order) : nx{_nx}, ny{_ny}, nz{_nx}, order(_order) {};
    void setnpx(int _npx, int _npy, int _npz)
    {
        npx = _npx;
        npy = _npy;
        npz = _npz;
    }
    void setnpx(int _npx)
    {
        npx = _npx;
        npy = _npx;
        npz = _npx;
    }
    void runPKernel(std::vector<std::vector<float>> &, std::map<std::string, std::vector<int>> &, std::vector<std::vector<float>> &);
    void createMemory(int &, int &, int &, float sigma, float Dq);

    std::vector<std::vector<float>> getSaxs();

    ~saxsKernel();

private:
    int size;
    int order;
    int npx, npy, npz;
    int nx, ny, nz, nnx, nny, nnz;
    int numParticles;
    float sigma;
    float bin_size;
    int num_bins;
    thrust::device_vector<float> d_moduleX;
    thrust::device_vector<float> d_moduleY;
    thrust::device_vector<float> d_moduleZ;
    thrust::device_vector<float> d_grid;
    thrust::device_vector<float> d_gridSup;
    thrust::device_vector<cuFloatComplex> d_gridSupC;
    thrust::device_vector<float> d_histogram;
    thrust::device_vector<float> d_nhist;
    float *d_grid_ptr{nullptr};
    float *d_gridSup_ptr{nullptr};
    cuFloatComplex *d_gridSupC_ptr{nullptr};
    // Do bspmod
    float *d_moduleX_ptr{nullptr};
    float *d_moduleY_ptr{nullptr};
    float *d_moduleZ_ptr{nullptr};
    float *d_histogram_ptr{nullptr};
    float *d_nhist_ptr{nullptr};
    std::function<int(int, double)> borderBins = [](int nx, double shell) -> int
    {
        return static_cast<int>(shell * nx / 2);
    };

    std::vector<long long> generateMultiples(long long limit);
    long long findClosestProduct(int n, double sigma);
    friend __global__ void calculate_histogram(cuFloatComplex *d_array, float *d_histogram, float *nhist, float *oc, int nx, int ny, int nz,
                                               float bin_size, int num_bins);

    friend __global__ void modulusKernel(cuFloatComplex *grid_q, float *modX, float *modY, float *modZ,
                                         int numParticles, int nnx, int nny, int nnz);

    friend __global__ void scatterKernel(cuFloatComplex *grid_q, cuFloatComplex *grid_oq, float *oc,
                                         float *Scatter, int nnx, int nny, int nnz);
    friend __global__ void rhoKernel(float *xa, float *grid, int order,
                                     int numParticles, int nx, int ny, int nz);
    friend __global__ void superDensityKernel(float *d_grid, float *d_gridSup, float myDens,
                                              int nx, int ny, int nz, int nnx, int nny, int nnz);
    friend __global__ void zeroDensityKernel(float *d_grid, int size);
    friend __global__ void zeroDensityKernel(cuFloatComplex *d_grid, size_t size);
    friend __global__ void paddingKernel(float *grid, int nx, int ny, int nz, int dx, int dy, int dz, float *Dens, int *count);
};

#endif