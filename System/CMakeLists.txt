#[[
 * @file CMakeLists.txt
 * @brief Defines the build configuration for the System library.
 *
 * This CMakeLists.txt file sets the minimum required CMake version, configures the
 * build type, specifies the source files to be compiled into the System library,
 * and defines a compile definition for the Python source directory.
]]
# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE Cell.cpp AtomCounter.cpp)
add_library(system STATIC ${SOURCE})
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")
message(STATUS "System Include Directories: ${PROJECT_SOURCE_DIR}/include")
