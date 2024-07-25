# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE Atoms.cpp Python3.cpp Cell.cpp)
add_library(system STATIC ${SOURCE})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/Utilities")
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")
message(STATUS "System Include Directories: ${PROJECT_SOURCE_DIR}/include")
