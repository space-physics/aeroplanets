INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libprotoncosmics.la
libprotoncosmics_la_SOURCES= protocos.cpp protocos.hpp
libprotoncosmics_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../phflux/libphflux.la $(FLIBS)



if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test/test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la  libprotoncosmics.la ../phflux/libphflux.la   $(FLIBS)
mstest_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif


