INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libprotimpact.la
libprotimpact_la_SOURCES= proton.cpp proton.hpp fortran/dgphiv.f fortran/expogk.f
libprotimpact_la_LDFLAGS=--no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../phflux/libphflux.la $(BLAS_LIBS) $(LIBS) $(FLIBS)

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test/test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libprotimpact.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../phflux/libphflux.la  $(BLAS_LIBS) $(LIBS) $(FLIBS)
##mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@ 
TESTS= mstest

endif



