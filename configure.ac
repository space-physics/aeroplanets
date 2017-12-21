#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([AeroPlanets], [1.01], [guillaume.gronoff@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main/main.hpp])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4]) # Demande par libtools
# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_F77
AC_PROG_RANLIB
AC_F77_LIBRARY_LDFLAGS

#
AM_MAINTAINER_MODE
# Checks for libraries.


AX_BOOST_BASE([1.33])
AX_BOOST_SYSTEM()
AX_BOOST_PROGRAM_OPTIONS()
AX_BOOST_THREAD()
AX_BOOST_REGEX()
AX_BOOST_FILESYSTEM()


if ${HAVE_BOOST} ; then
		AC_MSG_NOTICE([Super pour boost])
else
		AC_MSG_ERROR([No boost library, please install first])
fi



if ${HAVE_BOOST_PROGRAM_OPTIONS} ; then
		AC_MSG_NOTICE([Super pour boost options])
else
		AC_MSG_ERROR([No boost options detected, please install first])
fi

if ${HAVE_BOOST_REGEX} ; then
		AC_MSG_NOTICE([Super pour boost regex ])
else
		AC_MSG_ERROR([No boost regex detected, please install first])
fi

if ${HAVE_BOOST_FILESYSTEM} ; then
	AC_MSG_NOTICE([Super pour boost filesystem ])
else
	AC_MSG_ERROR([No boost filesystem detected, please install first])
fi



if ${HAVE_BOOST_UNIT_TEST_FRAMEWORK} ; then
		AX_BOOST_UNIT_TEST_FRAMEWORK()
		AC_MSG_NOTICE([Super pour les tests])
		AM_CONDITIONAL(USE_TESTS, true)
else
		AC_MSG_NOTICE([No boost unit test detected, please install first if you want to use them!!!])
		AM_CONDITIONAL(USE_TESTS, false)
fi



AX_DEBUG()


AX_BLAS()


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
#AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_CHECK_TYPES([ptrdiff_t])


# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strstr])


AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/cppscixml/Makefile
                 src/phflux/Makefile
                 src/eflux/Makefile
                 src/leastsq/Makefile
                 src/minpack/Makefile
                 src/emit/Makefile
                 src/main/Makefile
                 src/chem/Makefile
                 src/math/Makefile
                 src/tinyxpath/Makefile
		 src/species/Makefile
		 src/photo/Makefile
                 src/proton/Makefile
                 src/protoncosmics/Makefile
		 src/elecimpact/Makefile
                 src/observ/Makefile
		 src/planet/Makefile
		 src/atmo/Makefile
		 utiles/Makefile
		 MonteCarlo/Makefile
                 src/paramfitting/Makefile
                 src/chemsimu/Makefile
                 src/computebranchphot/Makefile
	 	 doc/Makefile])
AC_OUTPUT
