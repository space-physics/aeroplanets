add_library(mathfun OBJECT chapman.cpp mathflux.cpp mathrandom.cpp umathfunction.cpp logging.cpp mathfunction.cpp mathstring.cpp)
target_link_libraries(mathfun PRIVATE Boost::regex)
target_include_directories(mathfun PRIVATE ${PROJECT_SOURCE_DIR}/src
  ${Boost_INCLUDE_DIR})  # for Mac