INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libchem.la
libchem_la_SOURCES= chem.cpp chem.hpp documentation.hpp reaclist.cpp reaclist.hpp reaction.cpp reaction.hpp calcdens.cpp calcdens.hpp

libchem_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 


if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libchem.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
#mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system @BOOST_SYSTEM_LIB@ 
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif



