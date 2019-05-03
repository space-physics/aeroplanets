add_library(species crs.cpp eleccrs.cpp path.cpp protoncrs.cpp shirai.cpp species.cpp)
target_include_directories(species PRIVATE ${CMAKE_SOURCE_DIR}/src)
