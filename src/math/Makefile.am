INCLUDES=-I@top_srcdir@
lib_LTLIBRARIES=libmathstring.la
libmathstring_la_SOURCES=mathstring.cpp mathstring.hpp mathfunction.hpp mathfunction.cpp  mathflux.cpp mathflux.hpp constphys.hpp mathrandom.cpp mathrandom.hpp logging.cpp logging.hpp chapman.cpp chapman.hpp umathfunction.cpp umathfunction.hpp
libmathstring_la_LDFLAGS=--no-undefined @BOOST_LDFLAGS@ @BOOST_REGEX_LIB@

if USE_TESTS
check_PROGRAMS=mstest
mstest_SOURCES=test/test.cpp
mstest_LDADD = @BOOST_UNIT_TEST_FRAMEWORK_LIB@  @BOOST_REGEX_LIB@ libmathstring.la
mstest_LDFLAGS = @BOOST_LDFLAGS@
TESTS= mstest

endif
EXTRA_DIST=./test/test.cpp
