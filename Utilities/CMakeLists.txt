#[[
 * @file CMakeLists.txt
 * @brief Defines the build configuration for the 'utils' static library.
 *
 * This CMakeLists.txt file sets the minimum required CMake version, configures the
 * build type, and defines the source files to be compiled into the 'utils' static
 * library.
]]
cmake_minimum_required( VERSION 3.18) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE BSpline.cpp BSpmod.cpp Scattering.cpp )
add_library(utils STATIC ${SOURCE})
