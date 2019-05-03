# doc/html/atmo_8hpp.html#details
add_library(atmo atmo.cpp neutralatmo.cpp)
target_link_libraries(atmo PRIVATE species chem photo proton emit elecimpact protocos planet scixml)
target_include_directories(atmo PRIVATE ${CMAKE_SOURCE_DIR}/src)