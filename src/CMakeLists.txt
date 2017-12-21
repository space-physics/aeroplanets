cmake_minimum_required (VERSION 2.8.12)
project(trans++ CXX Fortran)

include_directories(src/species/)

add_library(scixml cppscixml/scixml.cpp)

add_library(species species/species.cpp)
target_link_libraries(species scixml)

# doc/html/atmo_8hpp.html#details
add_executable(aero atmo/atmo.cpp)
target_link_libraries(aero species)
