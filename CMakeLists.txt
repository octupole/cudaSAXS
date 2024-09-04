#[[
 * @file CMakeLists.txt
 * @brief CMake configuration file for the cudaSAXS project.
 *
 * This CMake file sets up the build environment for the cudaSAXS project, which
 * includes the following:
 *
 * - Specifies the minimum required CMake version and the project name.
 * - Configures the C++ and CUDA standards to be used.
 * - Finds and includes the necessary dependencies, such as Python3, pybind11,
 *   CUDA, cuBLAS, and OpenMP.
 * - Includes the necessary directories for CUDA, Python3, and pybind11.
 * - Sets the CUDA architectures to be supported.
 * - Enables optimization flags for CUDA and C++ compilation.
 * - Adds subdirectories for the project's components (fmt, Utilities, System,
 *   Saxs, Exec).
 * - Specifies the target executable (cudaSAXS) and links the necessary
 *   libraries.
 * - Adds a compile definition for the Python source directory.
 * - Prints debugging messages for the include and library directories.
]]
cmake_minimum_required(VERSION 3.18)
#list(APPEND CMAKE_PREFIX_PATH "/opt/anaconda3/envs/openSAXS")
project(cudaSAXS LANGUAGES CXX C CUDA)
# Specify the CUDA and C++ standard
set(CMAKE_CXX_STANDARD 14)


# Find Python and pybind11

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/opt/anaconda3/envs/cudasaxs/share/cmake/pybind11")
find_package(pybind11 REQUIRED)


# Find CUDA toolkit
find_package(CUDAToolkit REQUIRED)
# Add cuBLAS library
find_library(CUBLAS_LIBRARIES cublas HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64 ${CUDA_TOOLKIT_ROOT_DIR}/lib)

# Find OpenMP
find_package(OpenMP REQUIRED)

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

# Include directories for CUDA
include_directories(${CUDA_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/Exec")
# Include directories for Python3 and Pybind11
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# Add CUDA library path manually if necessary
link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)

# Set the CUDA architectures you want to support
set(CMAKE_CUDA_ARCHITECTURES 60 70 75 80 86 89)  # For Pascal Volta, Turing, Ampere, Ampere and Ada Lovelace)

# Enable optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -diag-suppress 128 -O3 -rdc=true  -Xcompiler -O3")

add_subdirectory(fmt)
add_subdirectory(Utilities)
add_subdirectory(System)
add_subdirectory(Saxs)
add_subdirectory(Exec)

# Specify the target executable
add_executable(cudaSAXS cudaSAXS.cu )

target_include_directories(cudaSAXS PUBLIC "${PROJECT_SOURCE_DIR}/Exec ${PROJECT_SOURCE_DIR}/System ")

# Enable C++ and CUDA compilation
set_target_properties(cudaSAXS PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link libraries
target_link_libraries(cudaSAXS exec saxs system utils ${CUDA_LIBRARIES} cufft cublas fmt::fmt Python3::Python pybind11::module OpenMP::OpenMP_CXX)
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")

# Debugging messages
message(STATUS "Python3 Include Directories: ${Python3_INCLUDE_DIRS}")
message(STATUS "pybind11 Include Directories: ${pybind11_INCLUDE_DIRS}")
message(STATUS "CUDA Include Directories: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directories: ${CUDA_TOOLKIT_ROOT_DIR}/lib64")
