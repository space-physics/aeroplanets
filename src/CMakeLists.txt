add_subdirectory(math)

add_subdirectory(tinyxpath)

add_subdirectory(cppscixml)

add_subdirectory(elecimpact)

add_library(eflux OBJECT eflux/eflux.cpp)
target_include_directories(eflux PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(photo)
#-------------
add_subdirectory(observ)

add_subdirectory(emit)
#----------
add_library(phflux OBJECT phflux/phflux.cpp)
target_include_directories(phflux PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(dgphiv OBJECT proton/fortran/dgphiv.f proton/fortran/expogk.f)

add_library(proton OBJECT proton/proton.cpp)
target_link_libraries(proton PRIVATE phflux dgphiv ${BLAS_LIBRARIES})
target_include_directories(proton PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(planet)

#--- indp

add_subdirectory(species)

add_subdirectory(chem)

add_library(protocos OBJECT protoncosmics/protocos.cpp)
target_include_directories(protocos PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(atmo)
