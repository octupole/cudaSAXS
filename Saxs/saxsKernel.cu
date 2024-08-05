#include "saxsKernel.h"
#include "BSpmod.h"
#include "Scattering.h"
#include "opsfact.h"
#include <cuda_runtime.h> // Include CUDA runtime header
#include <cuComplex.h>
// Kernel to calculate |K| values and populate the histogram
__global__ void calculate_histogram(cuFloatComplex *d_array, float *d_histogram, float *d_nhist, float *oc, int nx, int ny, int nz,
                                    float bin_size, float qcut, int num_bins)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int npz = nz / 2 + 1;
    if (i < nx && j < ny && k < npz)
    { // Only consider the upper half in z-direction

        int nfx = (nx % 2 == 0) ? nx / 2 : nx / 2 + 1;
        int nfy = (ny % 2 == 0) ? ny / 2 : ny / 2 + 1;
        int nfz = (nz % 2 == 0) ? nz / 2 : nz / 2 + 1;

        int ia = (i < nfx) ? i : i - nx;
        int ja = (j < nfy) ? j : j - ny;
        int ka = (k < nfz) ? k : k - nz;
        int ib = i == 0 ? 0 : nx - i;
        int jb = j == 0 ? 0 : ny - j;
        float mw1, mw2, mw3, mw;
        mw1 = oc[XX * DIM + XX] * ia + oc[XX * DIM + YY] * ja + oc[XX * DIM + ZZ] * ka;
        mw1 = 2.0 * M_PI * mw1;
        mw2 = oc[YY * DIM + XX] * ia + oc[YY * DIM + YY] * ja + oc[YY * DIM + ZZ] * ka;
        mw2 = 2.0 * M_PI * mw2;
        mw3 = oc[ZZ * DIM + XX] * ia + oc[ZZ * DIM + YY] * ja + oc[ZZ * DIM + ZZ] * ka;
        mw3 = 2.0 * M_PI * mw3;
        mw = sqrtf(mw1 * mw1 + mw2 * mw2 + mw3 * mw3);
        if (mw > qcut)
            return;
        int h0 = static_cast<int>(mw / bin_size);
        int h1 = h0 + 1;
        cuFloatComplex v0;
        if (h0 < num_bins)
        {
            int idx = k + j * npz + i * npz * ny;
            int idbx = k + jb * npz + ib * npz * ny;
            v0 = d_array[idx];
            if (k != 0 && k != npz - 1)
            {
                auto v1 = d_array[idbx];
                v0 = cuCaddf(v0, v1);
                v0 = cuCmulf(v0, make_cuFloatComplex(0.5f, 0.0f));
            }
            atomicAdd(&d_histogram[h0], cuCrealf(v0));
            atomicAdd(&d_nhist[h0], 1.0f);
            if (h0 != 0)
            {
                atomicAdd(&d_histogram[h1], cuCrealf(v0));
                atomicAdd(&d_nhist[h1], 1.0f);
            }
        }
    }
}

/**
 * @brief Applies a modulus calculation to a grid of complex values.
 *
 * This kernel function calculates the modulus of each complex value in the input grid
 * and stores the result in the output grid.
 *
 * @param grid_q The input grid of complex values.
 * @param modX The modulus values for the x-dimension.
 * @param modY The modulus values for the y-dimension.
 * @param modZ The modulus values for the z-dimension.
 * @param numParticles The number of particles.
 * @param nnx The number of grid points in the x-dimension.
 * @param nny The number of grid points in the y-dimension.
 * @param nnz The number of grid points in the z-dimension.
 */
__global__ void modulusKernel(cuFloatComplex *grid_q, float *modX, float *modY, float *modZ,
                              int numParticles, int nnx, int nny, int nnz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int nnpz = nnz / 2 + 1;
    if (i < nnx && j < nny && k < nnpz)
    {
        int idx = k + j * nnpz + i * nnpz * nny;
        float bsp_i = modX[i];
        float bsp_j = modX[j];
        float bsp_k = modX[k];
        float bsp_ijk = bsp_i * bsp_j * bsp_k / (float)numParticles;
        cuFloatComplex bsp = make_cuComplex(bsp_ijk, 0.0f);
        grid_q[idx] = cuCmulf(cuConjf(grid_q[idx]), grid_q[idx]);
        grid_q[idx] = cuCmulf(grid_q[idx], bsp);
    }
}
/**
 * @brief Performs scattering calculations on a grid of complex values.
 *
 * This kernel function calculates the scattering contribution for each grid point
 * based on the provided scattering factors and the grid of complex values.
 *
 * @param grid_q The input grid of complex values.
 * @param grid_oq The output grid of complex values.
 * @param oc The orientation coefficients.
 * @param Scatter The scattering factors.
 * @param nnx The number of grid points in the x-dimension.
 * @param nny The number of grid points in the y-dimension.
 * @param nnz The number of grid points in the z-dimension.
 */
__global__ void scatterKernel(cuFloatComplex *grid_q, cuFloatComplex *grid_oq, float *oc,
                              float *Scatter, int nnx, int nny, int nnz, float qcut)
{

    // if (idx >= nx0 * ny0 * (nz0 / 2 + 1))
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int nfx = (nnx % 2 == 0) ? nnx / 2 : nnx / 2 + 1;
    int nfy = (nny % 2 == 0) ? nny / 2 : nny / 2 + 1;
    int nfz = (nnz % 2 == 0) ? nnz / 2 : nnz / 2 + 1;
    int nnpz = nnz / 2 + 1;
    if (i < nnx && j < nny && k < nnpz)
    {
        int idx = k + j * nnpz + i * nnpz * nny;

        opsfact ff;
        ff.allocate_device(Scatter);
        int ia = (i < nfx) ? i : i - nnx;
        int ja = (j < nfy) ? j : j - nny;
        int ka = (k < nfz) ? k : k - nnz;
        float mw1, mw2, mw3, mw;
        mw1 = oc[XX * DIM + XX] * ia + oc[XX * DIM + YY] * ja + oc[XX * DIM + ZZ] * ka;
        mw2 = oc[YY * DIM + XX] * ia + oc[YY * DIM + YY] * ja + oc[YY * DIM + ZZ] * ka;
        mw3 = oc[ZZ * DIM + XX] * ia + oc[ZZ * DIM + YY] * ja + oc[ZZ * DIM + ZZ] * ka;
        mw1 = 2.0 * M_PI * mw1;
        mw2 = 2.0 * M_PI * mw2;
        mw3 = 2.0 * M_PI * mw3;
        mw = sqrt(mw1 * mw1 + mw2 * mw2 + mw3 * mw3);
        if (mw > qcut)
            return;
        cuFloatComplex fq = make_cuComplex(ff(mw), 0.0f);
        cuFloatComplex mult = cuCmulf(fq, grid_q[idx]);
        grid_oq[idx] = cuCaddf(grid_oq[idx], mult);
    }
}
__global__ void zeroDensityKernel(float *d_grid, size_t size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < (int)size)
    {
        d_grid[idx] = 0.0f;
    }
}
__global__ void zeroDensityKernel(cuFloatComplex *d_grid, size_t size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < (int)size)
    {
        d_grid[idx] = make_cuComplex(0.0f, 0.0f);
    }
}

/**
 * @brief Computes the density contribution of each particle to the grid.
 *
 * This kernel function calculates the density contribution of each particle to the grid
 * using B-spline interpolation. It iterates over the grid points within the support
 * of the particle and adds the contribution to the corresponding grid points.
 *
 * @param xa The array of particle coordinates.
 * @param grid The grid to store the density contributions.
 * @param order The order of the B-spline interpolation.
 * @param numParticles The number of particles.
 * @param nx The number of grid points in the x-dimension.
 * @param ny The number of grid points in the y-dimension.
 * @param nz The number of grid points in the z-dimension.
 */
__global__ void rhoKernel(float *xa, float *grid, int order, int numParticles, int nx, int ny, int nz)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < numParticles)
    {
        Splines bsplineX;
        Splines bsplineY;
        Splines bsplineZ;

        int nx0 = static_cast<int>(nx);
        int ny0 = static_cast<int>(ny);
        int nz0 = static_cast<int>(nz);
        float x1, y1, z1, r1, s1, t1, gx, gy, gz;
        int mx, my, mz;

        x1 = xa[idx * DIM + XX];
        y1 = xa[idx * DIM + YY];
        z1 = xa[idx * DIM + ZZ];
        r1 = static_cast<float>(nx0 * (x1 - rint(x1 - 0.5)));
        s1 = static_cast<float>(ny0 * (y1 - rint(y1 - 0.5)));
        t1 = static_cast<float>(nz0 * (z1 - rint(z1 - 0.5)));
        mx = static_cast<int>(r1);
        my = static_cast<int>(s1);
        mz = static_cast<int>(t1);

        gx = r1 - static_cast<float>(mx);
        gy = s1 - static_cast<float>(my);
        gz = t1 - static_cast<float>(mz);
        spline splX = bsplineX(gx);
        spline splY = bsplineX(gy);
        spline splZ = bsplineX(gz);
        int i0 = mx - order;

        for (auto o = 0; o < order; o++)
        {
            int i = i0 + (nx0 - ((i0 >= 0) ? nx0 : -nx0)) / 2;

            int j0 = my - order;
            for (auto p = 0; p < order; p++)
            {
                int j = j0 + (ny0 - ((j0 >= 0) ? ny0 : -ny0)) / 2;

                int k0 = mz - order;
                for (auto q = 0; q < order; q++)
                {
                    int k = k0 + (nz0 - ((k0 >= 0) ? nz0 : -nz0)) / 2;
                    float fact_o = splX.x[o];
                    float fact_p = fact_o * splY.x[p];
                    float fact_q = fact_p * splZ.x[q];
                    int ig = k + j * nz0 + i * nz0 * ny0;
                    atomicAdd(&grid[ig], fact_q);
                    k0++;
                }
                j0++;
            }
            i0++;
        }
    }
}
/**
 * @brief Kernel function to initialize a 3D grid with a given density value.
 *
 * This kernel function is used to initialize a 3D grid with a given density value. The grid is represented as a 1D array, and the kernel function calculates the 1D index from the 3D coordinates of each grid point.
 *
 * @param d_grid Pointer to the 1D array representing the 3D grid.
 * @param myDens The density value to be assigned to the grid.
 * @param nx The size of the grid in the x-dimension.
 * @param ny The size of the grid in the y-dimension.
 * @param nz The size of the grid in the z-dimension.
 * @param nnx The size of the super-sampled grid in the x-dimension.
 * @param nny The size of the super-sampled grid in the y-dimension.
 * @param nnz The size of the super-sampled grid in the z-dimension.
 */
__global__ void superDensityKernel(float *d_grid, float *d_gridSup, float myDens, int nx, int ny, int nz, int nnx, int nny, int nnz)
{
    float N1 = (float)nnx / (float)nx;
    float N2 = (float)nny / (float)ny;
    float N3 = (float)nnz / (float)nz;

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    float summ0 = -myDens * (N1 * N2 * N3 - 1.0) / (N1 * N2 * N3);
    if (x < nnx && y < nny && z < nnz)
    {
        int idx_s = z + y * nnz + x * nnz * nny;
        d_gridSup[idx_s] = myDens;
        if (x < nx && y < ny && z < nz)
        {
            int idx = z + y * nz + x * nz * ny;
            d_gridSup[idx_s] = d_grid[idx];
        }
        d_gridSup[idx_s] += summ0;
    }
}

/**
 * @brief Performs padding on a 3D grid, computing the average density and count of points on the border.
 *
 * This CUDA kernel function performs padding on a 3D grid, computing the average density and count of points on the border of the grid. The grid is represented as a 1D array, and the kernel function calculates the 1D index from the 3D coordinates of each grid point.
 *
 * @param grid Pointer to the 1D array representing the 3D grid of floating-point values.
 * @param nx The size of the grid in the x-dimension.
 * @param ny The size of the grid in the y-dimension.
 * @param nz The size of the grid in the z-dimension.
 * @param dx The padding size in the x-dimension.
 * @param dy The padding size in the y-dimension.
 * @param dz The padding size in the z-dimension.
 * @param Dens Pointer to a device-side float variable to store the total density of the border points.
 * @param count Pointer to a device-side integer variable to store the count of border points.
 */
__global__ void paddingKernel(float *grid, int nx, int ny, int nz, int dx, int dy, int dz, float *Dens, int *count)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    int mx = nx - dx;
    int my = ny - dy;
    int mz = nz - dz;
    if (x < nx && y < ny && z < nz)
    {
        int idx = z + y * nz + x * nz * ny;
        bool cond1 = (x > dx && x < mx) && (y > dy && y < my) && (z > dz && z < mz);
        if (!cond1)
        {
            atomicAdd(&count[0], 1);
            atomicAdd(&Dens[0], grid[idx]);
        }
    }
}

/**
 * Processes a set of particles and computes their contribution to the SAXS intensity.
 *
 * This function iterates over a set of particles, transforms their coordinates based on the orientation matrix,
 * and computes their contribution to the SAXS intensity. It then performs padding, supersampling, and Fourier
 * transform operations on the density grid to compute the final SAXS intensity.
 *
 * @param coords A vector of particle coordinates.
 * @param index_map A map of particle indices, where the keys are particle types and the values are vectors of indices.
 * @param oc The orientation matrix.
 */
void saxsKernel::runPKernel(int frame, float Time, std::vector<std::vector<float>> &coords, std::map<std::string, std::vector<int>> &index_map, std::vector<std::vector<float>> &oc)
{
    static bool firstTime = true;
    // Cudaevents

    // to compute average density on the border
    if (firstTime)
    {
        this->resetHistogramParameters(oc);
        this->createMemory();
        this->writeBanner();
    }

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    int mx = borderBins(nx, SHELL);
    int my = borderBins(ny, SHELL);
    int mz = borderBins(nz, SHELL);
    float mySigma = (float)Options::nx / (float)Options::nnx;

    thrust::host_vector<float> h_oc(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
        {
            h_oc[i * DIM + j] = mySigma * oc[i][j];
        }

    thrust::device_vector<float> d_oc = h_oc;
    float *d_oc_ptr = thrust::raw_pointer_cast(d_oc.data());

    dim3 blockDim(npx, npy, npz);
    dim3 gridDim((nnx + blockDim.x - 1) / blockDim.x,
                 (nny + blockDim.y - 1) / blockDim.y,
                 (nnz + blockDim.z - 1) / blockDim.z);
    dim3 gridDim0((nx + blockDim.x - 1) / blockDim.x,
                  (ny + blockDim.y - 1) / blockDim.y,
                  (nz + blockDim.z - 1) / blockDim.z);
    const int THREADS_PER_BLOCK = 256;
    int numBlocksGrid = (d_grid.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    int numBlocksGridSuperC = (d_gridSupC.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    int numBlocksGridSuperAcc = (d_gridSupAcc.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    int numBlocksGridSuper = (d_gridSup.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    // zeroes the Sup density grid
    zeroDensityKernel<<<numBlocksGridSuperAcc, THREADS_PER_BLOCK>>>(d_gridSupAcc_ptr, d_gridSupAcc.size());

    int totParticles = 0;
    std::string formatted_string = fmt::format("--> Frame: {:<7}  Time Step: {:.2f} fs", frame, Time);

    // Print the formatted string
    std::cout << formatted_string << std::endl;

    for (const auto &pair : index_map)
    {
        cufftHandle plan;
        cufftPlan3d(&plan, nnx, nny, nnz, CUFFT_R2C);

        thrust::host_vector<float> h_Dens = {0.0f};
        thrust::host_vector<int> h_count = {0};
        thrust::device_vector<float> d_Dens = h_Dens;
        thrust::device_vector<int> d_count = h_count;
        std::string type = pair.first;
        std::vector<int> value = pair.second;
        std::vector<std::vector<float>> Particles;

        std::transform(value.begin(), value.end(), std::back_inserter(Particles), [&coords](int i)
                       { return coords[i]; });

        this->numParticles = Particles.size();
        totParticles += this->numParticles;
        // Allocate and copy particles to the device
        thrust::host_vector<float> h_particles(numParticles * 3);
        for (int i = 0; i < numParticles; ++i)
        {
            h_particles[i * 3] = oc[XX][XX] * Particles[i][XX] + oc[XX][YY] * Particles[i][YY] + oc[XX][ZZ] * Particles[i][ZZ];
            h_particles[i * 3 + 1] = oc[YY][XX] * Particles[i][XX] + oc[YY][YY] * Particles[i][YY] + oc[YY][ZZ] * Particles[i][ZZ];
            h_particles[i * 3 + 2] = oc[ZZ][XX] * Particles[i][XX] + oc[ZZ][YY] * Particles[i][YY] + oc[ZZ][ZZ] * Particles[i][ZZ];
        }

        thrust::device_vector<float> d_particles = h_particles;
        thrust::host_vector<float> h_scatter = Scattering::getScattering(type);
        thrust::device_vector<float> d_scatter = h_scatter;

        float *d_particles_ptr = thrust::raw_pointer_cast(d_particles.data());
        float *d_scatter_ptr = thrust::raw_pointer_cast(d_scatter.data());

        int numBlocks = (numParticles + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

        //    Kernels launch for the rhoKernel

        zeroDensityKernel<<<numBlocksGrid, THREADS_PER_BLOCK>>>(d_grid_ptr, d_grid.size());

        rhoKernel<<<numBlocks, THREADS_PER_BLOCK>>>(d_particles_ptr, d_grid_ptr, order,
                                                    numParticles, nx, ny, nz);

        // Synchronize the device
        cudaDeviceSynchronize();
        paddingKernel<<<gridDim0, blockDim>>>(d_grid_ptr, nx, ny, nz, mx, my, mz,
                                              thrust::raw_pointer_cast(d_Dens.data()),
                                              thrust::raw_pointer_cast(d_count.data()));
        // Synchronize the device
        cudaDeviceSynchronize();

        h_Dens = d_Dens;
        h_count = d_count;
        float myDens = h_Dens[0] / (float)h_count[0];
        // zeroes the Sup density grid
        zeroDensityKernel<<<numBlocksGridSuperC, THREADS_PER_BLOCK>>>(d_gridSupC_ptr, d_gridSupC.size());
        cudaDeviceSynchronize();

        superDensityKernel<<<gridDim, blockDim>>>(d_grid_ptr, d_gridSup_ptr, myDens, nx, ny, nz, nnx, nny, nnz);

        // Synchronize the device
        cudaDeviceSynchronize();

        cufftExecR2C(plan, d_gridSup_ptr, d_gridSupC_ptr);
        cudaDeviceSynchronize();

        // Synchronize the device
        scatterKernel<<<gridDim, blockDim>>>(d_gridSupC_ptr, d_gridSupAcc_ptr, d_oc_ptr, d_scatter_ptr, nnx, nny, nnz, kcut);
        cudaDeviceSynchronize();
    }
    modulusKernel<<<gridDim, blockDim>>>(d_gridSupAcc_ptr, d_moduleX_ptr, d_moduleY_ptr, d_moduleZ_ptr, totParticles, nnx, nny, nnz);
    // // Synchronize the device
    cudaDeviceSynchronize();
    calculate_histogram<<<gridDim, blockDim>>>(d_gridSupAcc_ptr, d_histogram_ptr, d_nhist_ptr, d_oc_ptr, nnx, nny, nnz,
                                               bin_size, kcut, num_bins);

    cudaDeviceSynchronize();
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    // Calculate the elapsed time in milliseconds
    float gpuElapsedTime;
    cudaEventElapsedTime(&gpuElapsedTime, start, stop);
    // std::cout << "GPU Elapsed Time: " << gpuElapsedTime << " ms" << std::endl;

    // Destroy the events
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    firstTime = false;
}
std::vector<std::vector<float>> saxsKernel::getSaxs()
{

    std::vector<std::vector<float>> saxs;
    thrust::host_vector<float> h_histogram = d_histogram;
    thrust::host_vector<float> h_nhist = d_nhist;
    for (auto o{0}; o < h_histogram.size(); o++)
    {
        if (h_nhist[o] != 0.0f)
        {
            vector<float> val = {o * this->bin_size, h_histogram[o] / h_nhist[o]};
            saxs.push_back(val);
        }
    }
    return saxs;
}

/**
 * @brief Creates the necessary memory for the SAXS computation.
 *
 * This function sets up the memory buffers and allocates memory for the SAXS computation.
 * It calculates the optimal grid sizes (nnx, nny, nnz) based on the original grid sizes (nx, ny, nz)
 * and the given sigma value. It then creates the necessary host and device memory buffers for the
 * grid, super-grid, and module data.
 *
 * @param[in,out] nnx The optimal x-dimension of the super-grid.
 * @param[in,out] nny The optimal y-dimension of the super-grid.
 * @param[in,out] nnz The optimal z-dimension of the super-grid.
 * @param[in] sigma The sigma value used to calculate the optimal grid sizes.
 */
void saxsKernel::createMemory()
{
    size_t nnpz = nnz / 2 + 1;

    this->bin_size = Options::Dq;
    this->kcut = Options::Qcut;

    this->num_bins = static_cast<int>(kcut / bin_size) + 1;
    thrust::host_vector<float> h_histogram(num_bins, 0.0f);
    thrust::host_vector<float> h_nhist(num_bins, 0.0f);

    d_histogram = h_histogram;
    d_nhist = h_nhist;
    d_histogram_ptr = thrust::raw_pointer_cast(d_histogram.data());
    d_nhist_ptr = thrust::raw_pointer_cast(d_nhist.data());
    BSpline::BSpmod *bsp_modx = new BSpline::BSpmod(nnx, nny, nnz);

    thrust::host_vector<float> h_moduleX = bsp_modx->ModX();
    thrust::host_vector<float> h_moduleY = bsp_modx->ModY();
    thrust::host_vector<float> h_moduleZ = bsp_modx->ModZ();

    d_moduleX = h_moduleX;
    d_moduleY = h_moduleY;
    d_moduleZ = h_moduleZ;
    d_moduleX_ptr = thrust::raw_pointer_cast(d_moduleX.data());
    d_moduleY_ptr = thrust::raw_pointer_cast(d_moduleY.data());
    d_moduleZ_ptr = thrust::raw_pointer_cast(d_moduleZ.data());

    thrust::host_vector<float> h_grid(nx * ny * nz);
    thrust::host_vector<float> h_gridSup(nnx * nny * nnz);
    thrust::host_vector<cuFloatComplex> h_gridSupC(nnx * nny * nnpz);
    thrust::host_vector<cuFloatComplex> h_gridSupAcc(nnx * nny * nnpz);

    d_grid = h_grid;
    d_gridSup = h_gridSup;
    d_gridSupC = h_gridSupC;
    d_gridSupAcc = h_gridSupAcc;

    d_grid_ptr = thrust::raw_pointer_cast(d_grid.data());
    d_gridSup_ptr = thrust::raw_pointer_cast(d_gridSup.data());
    d_gridSupC_ptr = thrust::raw_pointer_cast(d_gridSupC.data());
    d_gridSupAcc_ptr = thrust::raw_pointer_cast(d_gridSupAcc.data());
    // Do bspmod
}
/**
 * Generates a vector of multiples of 2, 3, 5, and 7 up to a given limit.
 *
 * This function generates all possible multiples of 2, 3, 5, and 7 up to the
 * specified limit, and returns them as a sorted, unique vector.
 *
 * @param limit The maximum value to generate multiples up to.
 * @return A vector of all multiples of 2, 3, 5, and 7 up to the given limit.
 */
// Function to generate multiples of 2, 3, 5, and 7 up to a given limit
std::vector<long long> saxsKernel::generateMultiples(long long limit)
{
    std::vector<long long> multiples;
    for (int a = 0; std::pow(2, a) <= limit; ++a)
    {
        for (int b = 0; std::pow(2, a) * std::pow(3, b) <= limit; ++b)
        {
            for (int c = 0; std::pow(2, a) * std::pow(3, b) * std::pow(5, c) <= limit; ++c)
            {
                for (int d = 0; std::pow(2, a) * std::pow(3, b) * std::pow(5, c) * std::pow(7, d) <= limit; ++d)
                {
                    long long multiple = std::pow(2, a) * std::pow(3, b) * std::pow(5, c) * std::pow(7, d);
                    if (multiple <= limit)
                    {
                        multiples.push_back(multiple);
                    }
                }
            }
        }
    }
    std::sort(multiples.begin(), multiples.end());
    multiples.erase(std::unique(multiples.begin(), multiples.end()), multiples.end());
    return multiples;
}

/**
 * Finds the closest integer to N * sigma that is obtainable by multiplying only 2, 3, 5, and 7.
 *
 * This function takes a target value N and a standard deviation sigma, and finds the closest integer
 * to N * sigma that can be expressed as a product of only the prime factors 2, 3, 5, and 7.
 *
 * @param n The target value N.
 * @param sigma The standard deviation.
 * @return The closest integer to N * sigma that is obtainable by multiplying only 2, 3, 5, and 7.
 */
// Function to find the closest integer to N * sigma that is obtainable by multiplying only 2, 3, 5, and 7
long long saxsKernel::findClosestProduct(int n, float sigma)
{
    long long target = std::round(n * sigma);
    long long limit = target * 2; // A generous limit for generating multiples
    std::vector<long long> multiples = generateMultiples(limit);

    long long closest = target;
    long long minDifference = std::numeric_limits<long long>::max();

    for (long long multiple : multiples)
    {
        long long difference = std::abs(multiple - target);
        if (difference < minDifference)
        {
            minDifference = difference;
            closest = multiple;
        }
    }

    return closest;
}
void saxsKernel::scaledCell()
{
    sigma = Options::sigma;
    if (Options::nnx == 0)
    {
        nnx = this->nnx = static_cast<int>(findClosestProduct(nx, sigma));
        nny = this->nny = static_cast<int>(findClosestProduct(ny, sigma));
        nnz = this->nnz = static_cast<int>(findClosestProduct(nz, sigma));
        Options::nnx = nnx;
        Options::nny = nny;
        Options::nnz = nnz;
    }
    else
    {
        this->nnx = Options::nnx;
        this->nny = Options::nny;
        this->nnz = Options::nnz;
    }
}
void saxsKernel::resetHistogramParameters(std::vector<std::vector<float>> &oc)
{

    auto qcut = Options::Qcut;
    auto dq = Options::Dq;
    int nfx{(nnx % 2 == 0) ? nnx / 2 : nnx / 2 + 1};
    int nfy{(nny % 2 == 0) ? nny / 2 : nny / 2 + 1};
    int nfz{(nnz % 2 == 0) ? nnz / 2 : nnz / 2 + 1};
    float argx{2.0f * (float)M_PI * oc[XX][XX] / sigma};
    float argy{2.0f * (float)M_PI * oc[YY][YY] / sigma};
    float argz{2.0f * (float)M_PI * oc[ZZ][ZZ] / sigma};

    std::vector<float> fx{(float)nfx - 1, (float)nfy - 1, (float)nfz - 1};

    vector<float> mydq0 = {argx, argy, argz, dq};
    vector<float> mycut0 = {argx * fx[XX], argy * fx[YY], argz * fx[ZZ], qcut};

    dq = (*std::max_element(mydq0.begin(), mydq0.end()));
    qcut = *std::min_element(mycut0.begin(), mycut0.end());
    if (qcut != Options::Qcut)
    {
        std::string formatted_string = fmt::format("----- Qcut had to be reset to {:.2f} from  {:.2f} ----", qcut, Options::Qcut);
        std::cout << "\n--------------------------------------------------\n";
        std::cout << formatted_string << "\n";
        std::cout << "--------------------------------------------------\n\n";

        Options::Qcut = qcut;
    }
    if (dq != Options::Dq)
    {
        std::string formatted_string = fmt::format("----- Dq had to be reset to {:.3f} from  {:.3f} ----", dq, Options::Dq);
        std::cout << "\n--------------------------------------------------\n";
        std::cout << formatted_string << "\n";
        std::cout << "--------------------------------------------------\n\n";

        Options::Dq = dq;
    }
}
void saxsKernel::writeBanner()
{
    std::string banner = fmt::format(
        "*************************************************\n"
        "* {:^40}      *\n"
        "* {:<19} {:>4} * {:>4} * {:>4}        *\n"
        "* {:<19} {:>4} * {:>4} * {:>4}        *\n"
        "* {:<10} {:>4}      {:<10} {:>4}          *\n"
        "* {:<10} {:>4.3f}     {:<10} {:>2.f}      *\n"
        "*************************************************\n\n",
        "Running cudaSAXS", "Cell Grid", Options::nx, Options::ny, Options::nz,
        "Supercell Grid", Options::nnx, Options::nny, Options::nnz, "Order",
        Options::order, "Sigma", Options::sigma, "Bin Size", Options::Dq, "Q Cutoff ", Options::Qcut);

    std::cout << banner;
}

saxsKernel::~saxsKernel()
{
}
