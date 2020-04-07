add_library(observ OBJECT geoobs.cpp geopath.cpp geopoint.cpp)
target_include_directories(observ PRIVATE ${PROJECT_SOURCE_DIR}/src)
