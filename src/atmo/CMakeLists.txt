# doc/html/atmo_8hpp.html#details
add_library(atmo OBJECT atmo.cpp neutralatmo.cpp)
target_link_libraries(atmo PRIVATE species chem photo proton emit elecimpact protocos planet scixml Boost::regex Boost::headers)
target_include_directories(atmo
    PRIVATE ${PROJECT_SOURCE_DIR}/src
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})