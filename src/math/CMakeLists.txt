add_library(mathfun chapman.cpp mathflux.cpp mathrandom.cpp umathfunction.cpp logging.cpp mathfunction.cpp mathstring.cpp)
target_link_libraries(mathfun PRIVATE Boost::regex)
target_include_directories(mathfun PRIVATE ${CMAKE_SOURCE_DIR}/src)