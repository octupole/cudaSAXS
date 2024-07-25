cmake_minimum_required(VERSION 3.10)
project(openSAXS LANGUAGES CXX CUDA)
# Specify the CUDA and C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

# Find Python and pybind11

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

# Find CUDA toolkit
find_package(CUDA REQUIRED)

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

# Include directories for CUDA
include_directories(${CUDA_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/System")
# Include directories for Python3 and Pybind11
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# Add CUDA library path manually if necessary
link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)

# Set the CUDA architectures you want to support
set(CMAKE_CUDA_ARCHITECTURES 80)

# Enable optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ")
add_subdirectory(Cuda)
add_subdirectory(System)

# Specify the target executable
add_executable(openSAXS openSAXS.cpp)

target_include_directories(openSAXS PUBLIC "${PROJECT_SOURCE_DIR}/Cuda ${PROJECT_SOURCE_DIR}/System ")

# Enable C++ and CUDA compilation
set_target_properties(openSAXS PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link libraries
target_link_libraries(openSAXS cudafft3d system ${CUDA_LIBRARIES} cufft Python3::Python pybind11::module)
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")

# Debugging messages
message(STATUS "CUDA Include Directories: ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directories: ${CUDA_TOOLKIT_ROOT_DIR}/lib64")
