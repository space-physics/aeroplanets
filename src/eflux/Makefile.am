INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libeflux.la
libeflux_la_SOURCES=eflux.cpp eflux.hpp
libeflux_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la


if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test/test.cpp
mstest_LDADD =  @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libeflux.la ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif



EXTRA_DIST=./test/test.xml



