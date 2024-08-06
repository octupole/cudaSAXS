cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE BSpline.cpp BSpmod.cpp Scattering.cpp readTraj.cpp Fstream.cpp FstreamC.cpp)
add_library(utils STATIC ${SOURCE})
