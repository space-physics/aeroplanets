cmake_minimum_required (VERSION 3.13)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(trans++ CXX Fortran)

include(cmake/compilers.cmake)

find_package(BLAS REQUIRED)

if(WIN32 AND NOT BOOST_ROOT)
  # Boost 1.70 on Windows doesn't yet work with CMake 3.14
  set(BOOST_ROOT d:/boost_1_69_0)
  set(BOOST_INCLUDEDIR d:/boost_1_69_0/boost)
  set(BOOST_LIBRARYDIR d:/boost_1_69_0/lib64-msvc-14.1)
endif()
find_package(Boost REQUIRED COMPONENTS filesystem regex)

add_subdirectory(src)

add_executable(aero1d src/main/main.cpp)
target_link_libraries(aero1d atmo)
target_include_directories(aero1d PRIVATE src)

add_executable(aerochem src/chemsimu/main.cpp)
target_link_libraries(aerochem atmo)
target_include_directories(aerochem PRIVATE src)

add_executable(aerobranch src/computebranchphot/main.cpp)
target_link_libraries(aerobranch PRIVATE atmo)
target_include_directories(aerobranch PRIVATE src)