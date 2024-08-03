#include "saxsKernel.h"
#include "BSpmod.h"
#include "Scattering.h"
#include "opsfact.h"
#include <cuda_runtime.h> // Include CUDA runtime header
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
    int nx0 = static_cast<int>(nnx);
    int ny0 = static_cast<int>(nny);
    int nz0 = static_cast<int>(nnz);
    if (i < nx0 && j < ny0 && k < (nz0 / 2 + 1))
    {
        int idx = i + j * nx0 + k * nx0 * ny0;
        float bsp_i = modX[i];
        float bsp_j = modX[j];
        float bsp_k = modX[k];

        cuFloatComplex bsp = make_cuComplex(bsp_i * bsp_j * bsp_k, 0.0f);
        cuFloatComplex conj = cuConjf(grid_q[idx]);
        cuFloatComplex product = cuCmulf(conj, bsp);
        cuFloatComplex D = make_cuComplex(1.0f / (float)numParticles, 0.0f);
        grid_q[idx] = cuCmulf(cuCmulf(grid_q[idx], product), D);
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
                              float *Scatter, int nnx, int nny, int nnz)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int nx0 = static_cast<int>(nnx);
    int ny0 = static_cast<int>(nny);
    int nz0 = static_cast<int>(nnz);
    int nfx = (nx0 % 2 == 0) ? nx0 / 2 : nx0 / 2 + 1;
    int nfy = (ny0 % 2 == 0) ? ny0 / 2 : ny0 / 2 + 1;
    int nfz = (nz0 % 2 == 0) ? nz0 / 2 : nz0 / 2 + 1;
    if (i < nx0 && j < ny0 && k < (nz0 / 2 + 1))
    {
        int idx = i + j * nx0 + k * nx0 * ny0;
        // printf("i: %d, j: %d, k: %d\n", i, j, k);
        opsfact ff;
        ff.allocate_device(Scatter);
        int ia = (i < nfx) ? i : i - nx0;
        int ja = (j < nfy) ? j : j - ny0;
        int ka = (k < nfz) ? k : k - nz0;
        float mw1, mw2, mw3, mw;
        mw1 = oc[XX * DIM + XX] * ia + oc[XX * DIM + YY] * ja + oc[XX * DIM + ZZ] * ka;
        mw2 = oc[YY * DIM + XX] * ia + oc[YY * DIM + YY] * ja + oc[YY * DIM + ZZ] * ka;
        mw3 = oc[ZZ * DIM + XX] * ia + oc[ZZ * DIM + YY] * ja + oc[ZZ * DIM + ZZ] * ka;
        mw1 = 2.0 * M_PI * mw1;
        mw2 = 2.0 * M_PI * mw2;
        mw3 = 2.0 * M_PI * mw3;
        mw = sqrt(mw1 * mw1 + mw2 * mw2 + mw3 * mw3);
        cuFloatComplex fq = make_cuComplex(ff(mw), 0.0f);

        grid_oq[idx] = cuCaddf(grid_oq[idx], cuCmulf(fq, grid_q[idx]));
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

        float g_x[DIM][MAX_ORDER], g_dx[DIM][MAX_ORDER];
        Splines bsplineX;
        Splines bsplineY;
        Splines bsplineZ;
        bsplineX.allocate_device(&g_x[XX][0], &g_dx[XX][0]);
        bsplineY.allocate_device(&g_x[YY][0], &g_dx[YY][0]);
        bsplineZ.allocate_device(&g_x[ZZ][0], &g_dx[ZZ][0]);

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
                    int ig = i + j * nx0 + k * nx0 * ny0;
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

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    int nx0 = static_cast<int>(nnx);
    int ny0 = static_cast<int>(nny);
    int nz0 = static_cast<int>(nnz);
    if (x < nx0 && y < ny0 && z < nz0)
    {
        int idx_s = x + y * nx0 + z * nx0 * ny0;
        d_gridSup[idx_s] = myDens;
        if (x < nx && y < ny && z < nz)
        {
            int idx = x + y * nx + z * nx * ny;
            d_gridSup[idx_s] = d_grid[idx];
        }
    }
}
/**
 * @brief Initializes a 3D grid of complex numbers with zero values.
 *
 * This CUDA kernel function initializes a 3D grid of cuFloatComplex values to zero. The grid is represented as a 1D array, and the kernel function calculates the 1D index from the 3D coordinates of each grid point.
 *
 * @param d_grid Pointer to the 1D array representing the 3D grid of cuFloatComplex values.
 * @param nx The size of the grid in the x-dimension.
 * @param ny The size of the grid in the y-dimension.
 * @param nz The size of the grid in the z-dimension.
 */
__global__ void initializeIqKernel(cuFloatComplex *d_grid, int nx, int ny, int nz)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    int nx0 = static_cast<int>(nx);
    int ny0 = static_cast<int>(ny);
    int nz0 = static_cast<int>(nz);
    if (x < nx0 && y < ny0 && z < nz0)
    {
        int idx = x + y * nx0 + z * nx0 * ny0;
        d_grid[idx] = make_cuComplex(0.0f, 0.0f);
    }
}
/**
 * @brief Initializes a 3D grid of floating-point values to zero.
 *
 * This CUDA kernel function initializes a 3D grid of floating-point values to zero. The grid is represented as a 1D array, and the kernel function calculates the 1D index from the 3D coordinates of each grid point.
 *
 * @param d_grid Pointer to the 1D array representing the 3D grid of floating-point values.
 * @param nx The size of the grid in the x-dimension.
 * @param ny The size of the grid in the y-dimension.
 * @param nz The size of the grid in the z-dimension.
 */
__global__ void initializeDensityKernel(float *d_grid, int nx, int ny, int nz)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    int nx0 = static_cast<int>(nx);
    int ny0 = static_cast<int>(ny);
    int nz0 = static_cast<int>(nz);
    if (x < nx0 && y < ny0 && z < nz0)
    {
        int idx = x + y * nx0 + z * nx0 * ny0;
        d_grid[idx] = 0.0f;
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
    int nx0 = static_cast<int>(nx);
    int ny0 = static_cast<int>(ny);
    int nz0 = static_cast<int>(nz);
    int mx = nx - dx;
    int my = ny - dy;
    int mz = nz - dz;
    if (x < nx0 && y < ny0 && z < nz0)
    {
        int idx = x + y * nx0 + z * nx0 * ny0;
        bool cond1 = (x > dx && x < mx) && (y > dy && y < my) && (z > dz && z < mz);
        if (!cond1)
        {
            atomicAdd(count, 1);
            atomicAdd(Dens, grid[idx]);
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
void saxsKernel::runPKernel(std::vector<std::vector<float>> &coords, std::map<std::string, std::vector<int>> &index_map, std::vector<std::vector<float>> &oc)
{
    cufftHandle plan;
    cufftPlan3d(&plan, nnx, nny, nnz, CUFFT_R2C);
    // Cudaevents

    // cudaEvent_t start, stop;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord(start);

    // to compute average density on the border
    thrust::host_vector<float> h_Dens = {0.0f};
    thrust::host_vector<int> h_count = {0};
    thrust::device_vector<float> d_Dens = h_Dens;
    thrust::device_vector<int> d_count = h_count;
    int mx = borderBins(nx, SHELL);
    int my = borderBins(ny, SHELL);
    int mz = borderBins(nz, SHELL);
    thrust::host_vector<float> h_oc(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            h_oc[i * DIM + j] = oc[i][j];

    thrust::device_vector<float> d_oc = h_oc;
    float *d_oc_ptr = thrust::raw_pointer_cast(d_oc.data());

    dim3 blockDim(npx, npy, npz);
    dim3 gridDim((nnx + blockDim.x - 1) / blockDim.x,
                 (nny + blockDim.y - 1) / blockDim.y,
                 (nnz + blockDim.z - 1) / blockDim.z);
    auto nnpz = nnz / 2 + 1;
    initializeIqKernel<<<gridDim, blockDim>>>(d_gridSupC_ptr, nnx, nny, nnpz);

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
    for (const auto &pair : index_map)
    {
        std::string type = pair.first;
        std::vector<int> value = pair.second;
        std::vector<std::vector<float>> Particles;

        std::transform(value.begin(), value.end(), std::back_inserter(Particles), [&coords](int i)
                       { return coords[i]; });

        this->numParticles = Particles.size();

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

        const int THREADS_PER_BLOCK = 256;
        int numBlocks = (numParticles + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

        initializeDensityKernel<<<gridDim, blockDim>>>(d_grid_ptr, nx, ny, nz);

        // Synchronize the device
        cudaDeviceSynchronize();

        rhoKernel<<<numBlocks, THREADS_PER_BLOCK>>>(d_particles_ptr, d_grid_ptr, order,
                                                    numParticles, nx, ny, nz);
        // Synchronize the device
        cudaDeviceSynchronize();

        paddingKernel<<<gridDim, blockDim>>>(d_grid_ptr, nx, ny, nz, mx, my, mz,
                                             thrust::raw_pointer_cast(d_Dens.data()),
                                             thrust::raw_pointer_cast(d_count.data()));
        // Synchronize the device
        cudaDeviceSynchronize();

        h_Dens = d_Dens;
        h_count = d_count;
        float myDens = h_Dens[0] / (float)h_count[0];
        superDensityKernel<<<gridDim, blockDim>>>(d_grid_ptr, d_gridSup_ptr, myDens, nx, ny, nz, nnx, nny, nnz);

        // Synchronize the device
        cudaDeviceSynchronize();

        cufftExecR2C(plan, d_gridSup_ptr, (cuFloatComplex *)d_gridSup_ptr);

        // Synchronize the device
        cudaDeviceSynchronize();

        scatterKernel<<<gridDim, blockDim>>>((cuFloatComplex *)d_gridSup_ptr, d_gridSupC_ptr, d_oc_ptr, d_scatter_ptr, nnx, nny, nnz);

        // Synchronize the device
        cudaDeviceSynchronize();
        modulusKernel<<<gridDim, blockDim>>>(d_gridSupC_ptr, d_moduleX_ptr, d_moduleY_ptr, d_moduleZ_ptr, numParticles, nnx, nny, nnz);

        // Synchronize the device
        cudaDeviceSynchronize();
    }
    // cudaDeviceSynchronize();
    // cudaEventRecord(stop);
    // cudaEventSynchronize(stop);

    // // Calculate the elapsed time in milliseconds
    // float gpuElapsedTime;
    // cudaEventElapsedTime(&gpuElapsedTime, start, stop);

    // // Destroy the events
    // cudaEventDestroy(start);
    // cudaEventDestroy(stop);
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
void saxsKernel::createMemory(int &nnx, int &nny, int &nnz, float sigma)
{
    this->sigma = sigma;
    if (nnx == 0)
    {
        nnx = this->nnx = static_cast<int>(findClosestProduct(nx, sigma));
        nny = this->nny = static_cast<int>(findClosestProduct(ny, sigma));
        nnz = this->nnz = static_cast<int>(findClosestProduct(nz, sigma));
    }
    else
    {
        this->nnx = nnx;
        this->nny = nny;
        this->nnz = nnz;
    }

    size_t nnpz = nnz / 2 + 1;

    BSpline::BSpmod *bsp_modx = new BSpline::BSpmod(nx, ny, nz);
    std::cout << "Cell with nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
    std::cout << "SuperCell with nnx: " << nnx << " nny: " << nny << " nnz: " << nnz << std::endl;

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
    thrust::host_vector<float> h_gridSup(2 * nnx * nny * nnpz);
    thrust::host_vector<cuFloatComplex> h_gridSupC(nnx * nny * nnpz);

    d_grid = h_grid;
    d_gridSup = h_gridSup;
    d_gridSupC = h_gridSupC;

    d_grid_ptr = thrust::raw_pointer_cast(d_grid.data());
    d_gridSup_ptr = thrust::raw_pointer_cast(d_gridSup.data());
    d_gridSupC_ptr = thrust::raw_pointer_cast(d_gridSupC.data());
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
long long saxsKernel::findClosestProduct(int n, double sigma)
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

saxsKernel::~saxsKernel()
{
}
