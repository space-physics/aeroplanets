cmake_minimum_required (VERSION 2.8.12)
project(trans++ CXX Fortran)

include_directories(.)

add_compile_options(-g)

find_package(BLAS REQUIRED)

find_package(Boost REQUIRED COMPONENTS system filesystem regex unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})
#----------------

file(GLOB mthsrc math/*.cpp)
add_library(mathfun ${mthsrc})
target_link_libraries(mathfun PRIVATE boost_regex)

file(GLOB txpsrc tinyxpath/*.cpp)
add_library(tinyxpath ${txpsrc})

file(GLOB xmlsrc cppscixml/*.cpp)
add_library(scixml ${xmlsrc})
target_link_libraries(scixml PRIVATE mathfun tinyxpath boost_filesystem boost_unit_test_framework)
#---------
add_library(disort elecimpact/fortran/disort.f elecimpact/fortran/linsub.f)

add_library(elecimpact elecimpact/electronionization.cpp)
target_link_libraries(elecimpact PRIVATE disort)
#-----------
add_library(eflux eflux/eflux.cpp)

file(GLOB phtsrc photo/*.cpp)
add_library(photo ${phtsrc})
target_link_libraries(photo PRIVATE eflux)
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
#----------
add_library(iri planet/Earth/iris12.f planet/Earth/irif12.f planet/Earth/msis90.f planet/Earth/cira86_pour_iri12.f)
add_library(earth planet/Earth/irimsis.cpp)
target_link_libraries(earth PRIVATE iri)

add_library(venus planet/Venus/venssubs_V05.f planet/Venus/pvatmos.f)

add_library(titan planet/Titan/titandata.cpp)

file(GLOB marsrc planet/Mars/*.cpp)
add_library(mars ${marsrc})

file(GLOB pltsrc planet/*.cpp)
add_library(planet ${pltsrc})
target_link_libraries(planet PRIVATE earth venus titan mars)

#--- indp

file(GLOB spcsrc species/*.cpp)
add_library(species ${spcsrc})

file(GLOB chmsrc chem/*.cpp)
add_library(chem ${chmsrc})

add_library(protocos protoncosmics/protocos.cpp)

# doc/html/atmo_8hpp.html#details
add_library(atmo atmo/atmo.cpp atmo/neutralatmo.cpp)
target_link_libraries(atmo PRIVATE species chem photo proton emit elecimpact protocos planet scixml boost_system)


#--- Executables.
add_executable(aero1d main/main.cpp)
target_link_libraries(aero1d atmo)

add_executable(aerochem chemsimu/main.cpp)
target_link_libraries(aerochem atmo)

add_executable(aerobranch computebranchphot/main.cpp)
target_link_libraries(aerobranch atmo)
