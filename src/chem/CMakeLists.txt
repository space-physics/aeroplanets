add_library(chem calcdens.cpp chem.cpp reaclist.cpp reaction.cpp)
target_include_directories(chem PRIVATE ${CMAKE_SOURCE_DIR}/src)
