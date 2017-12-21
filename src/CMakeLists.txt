cmake_minimum_required (VERSION 2.8.12)
project(trans++ CXX Fortran)

include_directories(.)

find_package(BLAS REQUIRED)

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
#---------
add_library(disort elecimpact/fortran/disort.f elecimpact/fortran/linsub.f)

add_library(elecimpact elecimpact/electronionization.cpp)
target_link_libraries(elecimpact disort)
#-----------
add_library(eflux eflux/eflux.cpp)

file(GLOB phtsrc photo/*.cpp)
add_library(photo ${phtsrc})
target_link_libraries(photo eflux)
#-------------
file(GLOB obssrc observ/*.cpp)
add_library(observ ${obssrc})

file(GLOB emtsrc emit/*.cpp)
add_library(emit ${emtsrc})
target_link_libraries(emit PRIVATE observ)
#----------
add_library(phflux phflux/phflux.cpp)

add_library(dgphiv proton/fortran/dgphiv.f proton/fortran/expogk.f)

add_library(proton proton/proton.cpp)
target_link_libraries(proton PRIVATE phflux dgphiv blas)
#--- indp

file(GLOB spcsrc species/*.cpp)
add_library(species ${spcsrc})

file(GLOB chmsrc chem/*.cpp)
add_library(chem ${chmsrc})


# doc/html/atmo_8hpp.html#details
add_executable(aero atmo/atmo.cpp atmo/neutralatmo.cpp)
target_link_libraries(aero species chem photo proton emit elecimpact scixml boost_system)
