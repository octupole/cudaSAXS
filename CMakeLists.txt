cmake_minimum_required(VERSION 3.5...3.20)
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
# Find OpenMP
find_package(OpenMP REQUIRED)

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

# Include directories for CUDA
include_directories(${CUDA_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/System")
include_directories("${PROJECT_SOURCE_DIR}/cudaKernels")
include_directories("${PROJECT_SOURCE_DIR}/Saxs")
include_directories("${PROJECT_SOURCE_DIR}/Utilities")
include_directories("${PROJECT_SOURCE_DIR}/Exec")
include_directories("${PROJECT_SOURCE_DIR}/fmt/include")
include_directories("${PROJECT_SOURCE_DIR}/xdr/include")
# Include directories for Python3 and Pybind11
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# Add CUDA library path manually if necessary
link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)

# Set the CUDA architectures you want to support
set(CMAKE_CUDA_ARCHITECTURES 89)
# Enable optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -diag-suppress 128 -O3 -rdc=true -use_fast_math --generate-code arch=compute_89,code=[sm_89,compute_89] -Xcompiler -O3")

add_subdirectory(fmt)
add_subdirectory(Utilities)
add_subdirectory(System)
add_subdirectory(cudaKernels)
add_subdirectory(Saxs)
add_subdirectory(xdr)
add_subdirectory(Exec)

# Specify the target executable
add_executable(cudaSAXS cudaSAXS.cu )

target_include_directories(cudaSAXS PUBLIC "${PROJECT_SOURCE_DIR}/Exec ${PROJECT_SOURCE_DIR}/System ${PROJECT_SOURCE_DIR}/cudaKernels ")

# Enable C++ and CUDA compilation
set_target_properties(cudaSAXS PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link libraries
target_link_libraries(cudaSAXS cudaKernels exec saxs system utils xdrfile ${CUDA_LIBRARIES} cufft fmt::fmt Python3::Python pybind11::module OpenMP::OpenMP_CXX)
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")

# Debugging messages
message(STATUS "Python3 Include Directories: ${Python3_INCLUDE_DIRS}")
message(STATUS "pybind11 Include Directories: ${pybind11_INCLUDE_DIRS}")
message(STATUS "CUDA Include Directories: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directories: ${CUDA_TOOLKIT_ROOT_DIR}/lib64")
