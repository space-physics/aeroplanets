INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
bin_PROGRAMS=aerochemsimu
aerochemsimu_SOURCES=main.hpp main.cpp
##aerochemsimu_LDFLAGS=@BOOST_LDFLAGS@ -lboost_system $(FLIBS)  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../photo/libphoto.la  ../planet/libplanet.la ../atmo/libatmo.la ../elecimpact/libelecimpact.la ../phflux/libphflux.la  ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la
aerochemsimu_LDFLAGS=@BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@ $(FLIBS)  ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../photo/libphoto.la  ../planet/libplanet.la ../atmo/libatmo.la ../elecimpact/libelecimpact.la ../phflux/libphflux.la  ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la
