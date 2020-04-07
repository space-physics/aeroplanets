add_library(emit OBJECT emission.cpp emit.cpp emitlist.cpp)
target_link_libraries(emit PRIVATE observ scixml)
target_include_directories(emit
  PRIVATE ${PROJECT_SOURCE_DIR}/src
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})