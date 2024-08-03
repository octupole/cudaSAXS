cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE Splines.cu saxsKernel.cu opsfact.cu)
add_library(saxs STATIC ${SOURCE})
