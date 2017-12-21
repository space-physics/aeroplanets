INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libatmo.la
libatmo_la_SOURCES=atmo.hpp atmo.cpp  neutralatmo.cpp neutralatmo.hpp 
libatmo_la_LDFLAGS=$(FLIBS) --no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../photo/libphoto.la  ../planet/libplanet.la ../eflux/libeflux.la ../phflux/libphflux.la ../elecimpact/libelecimpact.la ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
mstest_LDADD = $(FLIBS) @BOOST_UNIT_TEST_FRAMEWORK_LIB@ libatmo.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../photo/libphoto.la ../planet/libplanet.la ../eflux/libeflux.la ../phflux/libphflux.la ../elecimpact/libelecimpact.la  ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la
## mstest_LDFLAGS = @BOOST_LDFLAGS@ -lboost_system
mstest_LDFLAGS = @BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@
TESTS= mstest

endif


EXTRA_DIST=./README ./test.cpp ./test/test.xml ./test/init_species.xml ./test/crs/photoionisation.xml
