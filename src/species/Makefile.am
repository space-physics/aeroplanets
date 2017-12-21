INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libspecies.la
libspecies_la_SOURCES=path.hpp path.cpp species.hpp species.cpp crs.cpp crs.hpp documentation.hpp eleccrs.cpp eleccrs.hpp shirai.cpp shirai.hpp protoncrs.cpp protoncrs.hpp
libspecies_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la 
if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@   libspecies.la
##mstest_LDFLAGS = @BOOST_LDFLAGS@  -lboost_system ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la
mstest_LDFLAGS = @BOOST_LDFLAGS@   @BOOST_SYSTEM_LIB@ ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la
TESTS= mstest

endif




EXTRA_DIST=./README ./test.cpp ./crs_test.xml
