INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libemit.la
libemit_la_SOURCES= emission.cpp emission.hpp emit.cpp emit.hpp emitlist.cpp emitlist.hpp

libemit_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../planet/libplanet.la ../observ/libobserv.la ../chem/libchem.la



if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
##mstest_LDADD = -lboost_system @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libemit.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../planet/libplanet.la  ../observ/libobserv.la ../chem/libchem.la
mstest_LDADD =  @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libemit.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../planet/libplanet.la  ../observ/libobserv.la ../chem/libchem.la
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@
TESTS= mstest

endif








