add_library(observ geoobs.cpp geopath.cpp geopoint.cpp)
target_include_directories(observ PRIVATE ${CMAKE_SOURCE_DIR}/src)
