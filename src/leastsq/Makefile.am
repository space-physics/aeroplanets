INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=liblsq.la
liblsq_la_SOURCES= least.cpp least.hpp
liblsq_la_LDFLAGS= -no-undefined ../math/libmathstring.la ../minpack/libminpack.la

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=./test/test.cpp
mstest_LDADD=  @BOOST_UNIT_TEST_FRAMEWORK_LIB@ liblsq.la ../math/libmathstring.la ../minpack/libminpack.la
##mstest_LDFLAGS= -no-undefined @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS= -no-undefined @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@
TESTS=mstest
endif



