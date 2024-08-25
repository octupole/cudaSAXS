#[[
 * @file Exec/CMakeLists.txt
 * @brief Defines the build configuration for the 'exec' library.
 *
 * This CMakeLists.txt file sets the minimum required CMake version, configures the build type, defines the source files for the 'exec' library, adds compile definitions, and includes necessary directories for the library.
 *
 * The 'exec' library is built as a static library and includes the following source files:
 * - RunSaxs.cu
 * - Options.cpp
 *
 * The library includes directories from the following project directories:
 * - Saxs
 * - System
 * - pysrc
 * - fmt/include
 *
 * The PY_SOURCE_DIR compile definition is also set to the pysrc directory within the project source directory.
 ]]
# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE RunSaxs.cu Options.cpp)
add_library(exec STATIC ${SOURCE})
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")
target_include_directories(exec PUBLIC "${PROJECT_SOURCE_DIR}/Saxs" "${PROJECT_SOURCE_DIR}/System" "${PROJECT_SOURCE_DIR}/pysrc"  "${PROJECT_SOURCE_DIR}/fmt/include")

