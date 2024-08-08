cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE saxsKernel.cu Splines.cu opsfact.cu saxsDeviceKernels.cu)
add_library(saxs STATIC ${SOURCE})
target_include_directories(saxs PRIVATE "${PROJECT_SOURCE_DIR}/Exec" "${PROJECT_SOURCE_DIR}/fmt/include" "${PROJECT_SOURCE_DIR}/Utilities" "${PROJECT_SOURCE_DIR}/cudaKernels")
