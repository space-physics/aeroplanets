add_library(scixml scixml.cpp)
target_include_directories(scixml PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(scixml PRIVATE mathfun tinyxpath Boost::filesystem)
