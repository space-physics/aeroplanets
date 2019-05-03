add_subdirectory(math)

add_subdirectory(tinyxpath)

add_subdirectory(cppscixml)

add_subdirectory(elecimpact)

add_library(eflux eflux/eflux.cpp)
target_include_directories(eflux PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(photo)
#-------------
add_subdirectory(observ)

add_subdirectory(emit)
#----------
add_library(phflux phflux/phflux.cpp)
target_include_directories(phflux PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(dgphiv proton/fortran/dgphiv.f proton/fortran/expogk.f)

add_library(proton proton/proton.cpp)
target_link_libraries(proton PRIVATE phflux dgphiv ${BLAS_LIBRARIES})
target_include_directories(proton PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(planet)

#--- indp

add_subdirectory(species)

add_subdirectory(chem)

add_library(protocos protoncosmics/protocos.cpp)
target_include_directories(protocos PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(atmo)

#--- Executables.
add_executable(aero1d main/main.cpp)
target_link_libraries(aero1d atmo)
target_include_directories(aero1d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(aerochem chemsimu/main.cpp)
target_link_libraries(aerochem atmo)
target_include_directories(aerochem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(aerobranch computebranchphot/main.cpp)
target_link_libraries(aerobranch PRIVATE atmo)
target_include_directories(aerobranch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
