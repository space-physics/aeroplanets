add_library(species OBJECT crs.cpp eleccrs.cpp path.cpp protoncrs.cpp shirai.cpp species.cpp)
target_link_libraries(species PRIVATE scixml)
target_include_directories(species
    PRIVATE ${PROJECT_SOURCE_DIR}/src
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
