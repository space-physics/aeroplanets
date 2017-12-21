INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
lib_LTLIBRARIES=libplanet.la
libplanet_la_SOURCES= planet.cpp planet.hpp Venus/pvatmos.f Venus/venssubs_V05.f allplanets.hpp mars.cpp mars.hpp venus.cpp venus.hpp Mars/marsbinary.hpp Mars/marsbinary.cpp Mars/datamarinerviking.hpp Mars/datamarinerviking.cpp Mars/marsatmotim.cpp Mars/marsatmotim.hpp titan.cpp titan.hpp Titan/titandata.cpp Titan/titandata.hpp Earth/cira86_pour_iri12.f Earth/irif12.f Earth/irimsis.cpp Earth/irimsis.hpp Earth/iris12.f Earth/msis90.f earth.cpp earth.hpp anonym.cpp anonym.hpp


libplanet_la_LDFLAGS=$(FLIBS) --no-undefined ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 



if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test.cpp
##mstest_LDADD = $(FLIBS) @BOOST_UNIT_TEST_FRAMEWORK_LIB@  libplanet.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 
mstest_LDADD = $(FLIBS) @BOOST_UNIT_TEST_FRAMEWORK_LIB@  -lboost_system libplanet.la  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la 
mstest_LDFLAGS = @BOOST_LDFLAGS@
TESTS= mstest

endif



EXTRA_DIST= ./test.cpp ./test/test.xml
