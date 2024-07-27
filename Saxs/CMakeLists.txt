# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE RunSaxs.cpp)
add_library(saxs STATIC ${SOURCE})
add_compile_definitions(PY_SOURCE_DIR="${CMAKE_SOURCE_DIR}/pysrc")
target_include_directories(saxs PUBLIC "${PROJECT_SOURCE_DIR}/Utilities ${PROJECT_SOURCE_DIR}/System ${PROJECT_SOURCE_DIR}/include")