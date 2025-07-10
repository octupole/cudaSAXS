#[[
 * @file CMakeLists.txt
 * @brief Defines the build configuration for the SAXS library.
 *
 * This CMakeLists.txt file sets the minimum required CMake version, configures the build type, and defines the source files to be compiled into the SAXS static library.
 * The library is configured to include directories from the Exec, fmt, Utilities, and cudaKernels directories.
 ]]
cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE saxsKernel.cu Splines.cu opsfact.cu saxsDeviceKernels.cu)
add_library(saxs STATIC ${SOURCE})
target_include_directories(saxs BEFORE PRIVATE "${PROJECT_SOURCE_DIR}/Exec" "${PROJECT_SOURCE_DIR}/fmt/include" "${PROJECT_SOURCE_DIR}/Utilities" "${PROJECT_SOURCE_DIR}/cudaKernels")
