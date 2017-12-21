INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libscixml.la
libscixml_la_SOURCES=scixml.hpp scixml.cpp error.hpp
libscixml_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la @BOOST_FILESYSTEM_LIB@

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
##mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@  -lboost_system libscixml.la ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la @BOOST_FILESYSTEM_LIB@
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@  libscixml.la ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la @BOOST_FILESYSTEM_LIB@
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif
EXTRA_DIST=./test.cpp ./test.xml ./test2.xml ./test3.xml




