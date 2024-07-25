# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE RhoSaxsBSP.cpp RhoSaxs.cpp )
add_library(saxs STATIC ${SOURCE})
target_include_directories(saxs PUBLIC "${PROJECT_SOURCE_DIR}/Utilities"" ${PROJECT_SOURCE_DIR}/include")