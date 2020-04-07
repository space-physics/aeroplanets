add_library(scixml OBJECT scixml.cpp)
target_include_directories(scixml
    PRIVATE ${PROJECT_SOURCE_DIR}/src
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(scixml PRIVATE mathfun tinyxpath Boost::filesystem)
