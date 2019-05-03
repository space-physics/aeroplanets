cmake_minimum_required (VERSION 3.13)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(trans++ CXX Fortran)

find_package(BLAS REQUIRED)

find_package(Boost REQUIRED COMPONENTS filesystem regex)

add_subdirectory(src)