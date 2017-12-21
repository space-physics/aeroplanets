INCLUDES=-I@top_srcdir@ -I@top_srcdir@/src
bin_PROGRAMS=aero1d
aero1d_SOURCES=main.hpp main.cpp documentation.hpp timeline.hpp
##aero1d_LDFLAGS=@BOOST_LDFLAGS@ -lboost_system $(BLAS_LIBS) $(LIBS) $(FLIBS) ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../photo/libphoto.la  ../planet/libplanet.la ../atmo/libatmo.la ../elecimpact/libelecimpact.la ../phflux/libphflux.la  ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la
aero1d_LDFLAGS=@BOOST_LDFLAGS@  @BOOST_SYSTEM_LIB@ $(BLAS_LIBS) $(LIBS) $(FLIBS) ../math/libmathstring.la  ../tinyxpath/libtinyxpath.la ../cppscixml/libscixml.la ../species/libspecies.la ../eflux/libeflux.la ../photo/libphoto.la  ../planet/libplanet.la ../atmo/libatmo.la ../elecimpact/libelecimpact.la ../phflux/libphflux.la  ../proton/libprotimpact.la ../protoncosmics/libprotoncosmics.la  ../chem/libchem.la ../observ/libobserv.la ../emit/libemit.la






