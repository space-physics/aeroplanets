cmake_minimum_required (VERSION 3.0)
project(trans++ CXX Fortran)

add_compile_options(-g)

find_package(BLAS REQUIRED)

find_package(Boost REQUIRED COMPONENTS filesystem regex)
#----------------
add_library(mathfun math/chapman.cpp math/mathflux.cpp math/mathrandom.cpp math/umathfunction.cpp math/logging.cpp math/mathfunction.cpp math/mathstring.cpp)
target_link_libraries(mathfun PRIVATE Boost::regex)
target_include_directories(mathfun PRIVATE .)

add_library(tinyxpath action_store.cpp tinystr.cpp tinyxmlparser.cpp xpath_expression.cpp xpath_static.cpp lex_util.cpp tinyxml.cpp tokenlist.cpp xpath_processor.cpp xpath_stream.cpp node_set.cpp tinyxmlerror.cpp xml_util.cpp xpath_stack.cpp xpath_syntax.cpp)

add_library(scixml cppscixml/scixml.cpp cppscixml/test.cpp)
target_include_directories(scixml PRIVATE ./)
target_link_libraries(scixml PRIVATE mathfun tinyxpath Boost::filesystem)
#---------
add_library(disort elecimpact/fortran/disort.f elecimpact/fortran/linsub.f)
target_compile_options(disort PRIVATE -w)

add_library(elecimpact elecimpact/electronionization.cpp)
target_link_libraries(elecimpact PRIVATE disort)
target_include_directories(elecimpact PRIVATE ./)
#-----------
add_library(eflux eflux/eflux.cpp)
target_include_directories(eflux PRIVATE ./)

file(GLOB phtsrc photo/*.cpp)
add_library(photo ${phtsrc})
target_link_libraries(photo PRIVATE eflux)
target_include_directories(photo PRIVATE ./)
#-------------
file(GLOB obssrc observ/*.cpp)
add_library(observ ${obssrc})
target_include_directories(observ PRIVATE ./)

file(GLOB emtsrc emit/*.cpp)
add_library(emit ${emtsrc})
target_link_libraries(emit PRIVATE observ)
target_include_directories(emit PRIVATE ./)
#----------
add_library(phflux phflux/phflux.cpp)
target_include_directories(phflux PRIVATE ./)

add_library(dgphiv proton/fortran/dgphiv.f proton/fortran/expogk.f)

add_library(proton proton/proton.cpp)
target_link_libraries(proton PRIVATE phflux dgphiv ${BLAS_LIBRARIES})
target_include_directories(proton PRIVATE ./)
#----------
add_library(iri planet/Earth/iris12.f planet/Earth/irif12.f planet/Earth/msis90.f planet/Earth/cira86_pour_iri12.f)
target_compile_options(iri PRIVATE -w)

add_library(earth planet/Earth/irimsis.cpp)
target_link_libraries(earth PRIVATE iri)
target_include_directories(earth PRIVATE ./)

add_library(venus planet/Venus/venssubs_V05.f planet/Venus/pvatmos.f)

add_library(titan planet/Titan/titandata.cpp)
target_include_directories(titan PRIVATE ./)

add_library(mars planet/Mars/datamarinerviking.cpp planet/Mars/marsatmotim.cpp planet/Mars/marsbinary.cpp)
target_include_directories(mars PRIVATE ./)

add_library(planet planet/anonym.cpp planet/earth.cpp planet/mars.cpp planet/planet.cpp planet/test.cpp planet/titan.cpp planet/venus.cpp)
target_link_libraries(planet PRIVATE earth venus titan mars)
target_include_directories(planet PRIVATE ./)

#--- indp

file(GLOB spcsrc species/*.cpp)
add_library(species ${spcsrc})
target_include_directories(species PRIVATE ./)

file(GLOB chmsrc chem/*.cpp)
add_library(chem ${chmsrc})
target_include_directories(chem PRIVATE ./)

add_library(protocos protoncosmics/protocos.cpp)
target_include_directories(protocos PRIVATE ./)

# doc/html/atmo_8hpp.html#details
add_library(atmo atmo/atmo.cpp atmo/neutralatmo.cpp)
target_link_libraries(atmo PRIVATE species chem photo proton emit elecimpact protocos planet scixml)
target_include_directories(atmo PRIVATE ./)


#--- Executables.
add_executable(aero1d main/main.cpp)
target_link_libraries(aero1d atmo)
target_include_directories(aero1d PRIVATE ./)

add_executable(aerochem chemsimu/main.cpp)
target_link_libraries(aerochem atmo)
target_include_directories(aerochem PRIVATE ./)

add_executable(aerobranch computebranchphot/main.cpp)
target_link_libraries(aerobranch PRIVATE atmo)
target_include_directories(aerobranch PRIVATE ./)
