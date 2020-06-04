add_library(disort OBJECT fortran/disort.f fortran/linsub.f)
target_compile_options(disort PRIVATE -w ${gcc10opts})

add_library(elecimpact OBJECT electronionization.cpp)
target_link_libraries(elecimpact PRIVATE disort)
target_include_directories(elecimpact PRIVATE ${PROJECT_SOURCE_DIR}/src)
