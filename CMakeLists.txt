cmake_minimum_required (VERSION 3.15)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(trans++
LANGUAGES C CXX Fortran)
enable_testing()

include(cmake/compilers.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

find_package(LAPACK REQUIRED)

# often need to set BOOST_ROOT and similar in your .bashrc or ~/boost.ps1 etc.
find_package(Boost REQUIRED COMPONENTS filesystem regex)

# compiling in parallel takes several GByte of RAM, let's throttle build parallelism if low RAM system
set_property(GLOBAL PROPERTY JOB_POOLS one_jobs=1 two_jobs=2 four_jobs=4)
cmake_host_system_information(RESULT _memfree QUERY AVAILABLE_PHYSICAL_MEMORY)
if(_memfree LESS 1000)
  set(CMAKE_JOB_POOL_COMPILE one_jobs)
elseif(_memfree LESS 2000)
  set(CMAKE_JOB_POOL_COMPILE two_jobs)
elseif(_memfree LESS 4000)
  set(CMAKE_JOB_POOL_COMPILE four_jobs)
endif()

add_subdirectory(src)

set(_intlibs atmo chem dgphiv disort earth eflux elecimpact emit iri mars mathfun observ phflux photo planet protocos proton scixml species tinyxpath titan venus)

add_executable(aero1d src/main/main.cpp)
target_link_libraries(aero1d ${_intlibs} Boost::headers Boost::regex)
target_include_directories(aero1d PRIVATE src)

add_executable(aerochem src/chemsimu/main.cpp)
target_link_libraries(aerochem ${_intlibs} Boost::headers Boost::regex)
target_include_directories(aerochem PRIVATE src)

add_executable(aerobranch src/computebranchphot/main.cpp)
target_link_libraries(aerobranch PRIVATE ${_intlibs} Boost::headers Boost::regex)
target_include_directories(aerobranch PRIVATE src)

# --- tests

add_test(NAME Earth
  COMMAND $<TARGET_FILE:aero1d> ${PROJECT_SOURCE_DIR}/data/Earth/AuroraEarthFairbanks.xml)