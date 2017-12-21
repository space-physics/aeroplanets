INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libphoto.la
libphoto_la_SOURCES=flux_model.cpp flux_model.hpp photoionization.cpp photoionization.hpp documentation.hpp lines.cpp lines.hpp
libphoto_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../phflux/libphflux.la

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libphoto.la   ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../phflux/libphflux.la
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif




