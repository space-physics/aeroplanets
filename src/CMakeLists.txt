cmake_minimum_required (VERSION 2.8.12)
project(trans++ CXX Fortran)

include_directories(.)

find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIR})
#----------------

file(GLOB mthsrc math/*.cpp)
add_library(mathfun ${mthsrc})
target_link_libraries(mathfun boost_regex)

file(GLOB txpsrc tinyxpath/*.cpp)
add_library(tinyxpath ${txpsrc})

file(GLOB xmlsrc cppscixml/*.cpp)
add_library(scixml ${xmlsrc})
target_link_libraries(scixml mathfun tinyxpath boost_filesystem boost_unit_test_framework)

file(GLOB spcsrc species/*.cpp)
add_library(species ${spcsrc})
#target_link_libraries(species scixml)

# doc/html/atmo_8hpp.html#details
add_executable(aero atmo/atmo.cpp)
target_link_libraries(aero species scixml boost_system)
