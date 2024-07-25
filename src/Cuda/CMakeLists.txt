cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE fft3d.cu)
add_library(cudafft3d STATIC ${SOURCE})
target_include_directories(cudafft3d PUBLIC "${PROJECT_SOURCE_DIR}/include")