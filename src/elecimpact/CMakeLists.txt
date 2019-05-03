add_library(disort fortran/disort.f fortran/linsub.f)
target_compile_options(disort PRIVATE -w)

add_library(elecimpact electronionization.cpp)
target_link_libraries(elecimpact PRIVATE disort)
target_include_directories(elecimpact PRIVATE ${CMAKE_SOURCE_DIR}/src)