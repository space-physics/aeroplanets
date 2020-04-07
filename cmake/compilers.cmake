set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)

#include(CheckIPOSupported)
#check_ipo_supported(RESULT lto_supported)
#if(lto_supported)
#  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
#endif()


if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
    add_compile_options(/arch:native)
    # /heap-arrays is necessary to avoid runtime errors in programs using this library
    string(APPEND CMAKE_Fortran_FLAGS " /heap-arrays")
  else()
    add_compile_options(-march=native)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native)
endif()
