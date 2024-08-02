cmake_minimum_required(VERSION 3.5...3.20)
#list(APPEND CMAKE_PREFIX_PATH "/opt/anaconda3/envs/openSAXS")
project(cudaSAXS LANGUAGES CXX CUDA)
# Specify the CUDA and C++ standard
set(CMAKE_CXX_STANDARD 14)


# Find Python and pybind11

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)


# Find CUDA toolkit
find_package(CUDAToolkit REQUIRED)

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

# Include directories for CUDA
include_directories(${CUDA_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/System")
include_directories("${PROJECT_SOURCE_DIR}/Saxs")
include_directories("${PROJECT_SOURCE_DIR}/Utilities")
# Include directories for Python3 and Pybind11
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# Add CUDA library path manually if necessary
link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)

# Set the CUDA architectures you want to support
set(CMAKE_CUDA_ARCHITECTURES 75)
# Enable optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -rdc=true --generate-code arch=compute_75,code=[sm_75,compute_75] -Xcompiler -O3")
add_subdirectory(Cuda)
add_subdirectory(System)
add_subdirectory(Saxs)

# Specify the target executable
add_executable(cudaSAXS cudaSAXS.cpp )

target_include_directories(cudaSAXS PUBLIC "${PROJECT_SOURCE_DIR}/Cuda ${PROJECT_SOURCE_DIR}/System ")

# Enable C++ and CUDA compilation
set_target_properties(cudaSAXS PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link libraries
target_link_libraries(cudaSAXS cudaLibs saxs system ${CUDA_LIBRARIES} cufft Python3::Python pybind11::module)
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")

# Debugging messages
message(STATUS "CUDA Include Directories: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directories: ${CUDA_TOOLKIT_ROOT_DIR}/lib64")
