add_library(photo OBJECT flux_model.cpp lines.cpp photoionization.cpp)
target_link_libraries(photo PRIVATE eflux scixml)
target_include_directories(photo PRIVATE ${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
