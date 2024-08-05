cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE saxsDeviceKernels.cu)
add_library(cudaKernels STATIC ${SOURCE})
