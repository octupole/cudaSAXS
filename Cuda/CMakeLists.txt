cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE CudaSaxsManager.cu)
add_library(cudaLibs STATIC ${SOURCE})
#target_include_directories(cudaLibs PUBLIC "${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Saxs ${PROJECT_SOURCE_DIR}/System ${PROJECT_SOURCE_DIR}/Utilities")