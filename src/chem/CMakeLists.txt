add_library(chem OBJECT calcdens.cpp chem.cpp reaclist.cpp reaction.cpp)
target_link_libraries(chem PRIVATE scixml)
target_include_directories(chem PRIVATE ${PROJECT_SOURCE_DIR}/src)
