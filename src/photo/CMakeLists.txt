add_library(photo flux_model.cpp lines.cpp photoionization.cpp test.cpp)
target_link_libraries(photo PRIVATE eflux)
target_include_directories(photo PRIVATE ${CMAKE_SOURCE_DIR}/src)