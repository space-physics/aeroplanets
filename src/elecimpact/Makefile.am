INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libelecimpact.la
libelecimpact_la_SOURCES= electronionization.cpp electronionization.hpp fortran/disort.f fortran/linsub.f
libelecimpact_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la ../phflux/libphflux.la $(FLIBS)

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test/test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libelecimpact.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la ../phflux/libphflux.la  $(FLIBS)
## mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@ @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif




