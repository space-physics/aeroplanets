set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)

#include(CheckIPOSupported)
#check_ipo_supported(RESULT lto_supported)
#if(lto_supported)
#  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
#endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  add_compile_options(
    $<$<IF:$<BOOL:${WIN32}>,/QxHost,-xHost>
    $<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native)
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  set(gcc10opts -fallow-argument-mismatch)
endif()
