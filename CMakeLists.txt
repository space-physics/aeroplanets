cmake_minimum_required (VERSION 3.15)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(trans++
LANGUAGES C CXX Fortran)

include(cmake/compilers.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

find_package(LAPACK REQUIRED)

# often need to set BOOST_ROOT and similar in your .bashrc or ~/boost.ps1 etc.
find_package(Boost REQUIRED COMPONENTS filesystem regex)

add_subdirectory(src)

set(_intlibs atmo chem dgphiv disort earth eflux elecimpact emit iri mars mathfun observ phflux photo planet protocos proton scixml species tinyxpath titan venus)

add_executable(aero1d src/main/main.cpp)
target_link_libraries(aero1d ${_intlibs})
target_include_directories(aero1d PRIVATE src)

add_executable(aerochem src/chemsimu/main.cpp)
target_link_libraries(aerochem ${_intlibs})
target_include_directories(aerochem PRIVATE src)

add_executable(aerobranch src/computebranchphot/main.cpp)
target_link_libraries(aerobranch PRIVATE ${_intlibs})
target_include_directories(aerobranch PRIVATE src)