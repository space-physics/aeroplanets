add_library(emit emission.cpp emit.cpp emitlist.cpp)
target_link_libraries(emit PRIVATE observ)
target_include_directories(emit PRIVATE ${CMAKE_SOURCE_DIR}/src)