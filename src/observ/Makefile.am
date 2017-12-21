INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libobserv.la
libobserv_la_SOURCES= geoobs.cpp geoobs.hpp geopath.cpp geopath.hpp geopoint.cpp geopoint.hpp documentation.hpp

libobserv_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../planet/libplanet.la


if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libobserv.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../planet/libplanet.la
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif




