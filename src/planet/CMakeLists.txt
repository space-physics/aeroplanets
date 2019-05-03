add_library(iri Earth/iris12.f Earth/irif12.f Earth/msis90.f Earth/cira86_pour_iri12.f)
target_compile_options(iri PRIVATE -w)

add_library(earth Earth/irimsis.cpp)
target_link_libraries(earth PRIVATE iri)
target_include_directories(earth PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_library(venus Venus/venssubs_V05.f Venus/pvatmos.f)

add_library(titan Titan/titandata.cpp)
target_include_directories(titan PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_library(mars Mars/datamarinerviking.cpp Mars/marsatmotim.cpp Mars/marsbinary.cpp)
target_include_directories(mars PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_library(planet anonym.cpp earth.cpp mars.cpp planet.cpp titan.cpp venus.cpp)
target_link_libraries(planet PRIVATE earth venus titan mars)
target_include_directories(planet PRIVATE ${CMAKE_SOURCE_DIR}/src)
