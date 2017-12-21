dnl Choppe sur le cours en gpl, puis pas mal modifie

AC_DEFUN([AX_DEBUG],
[
AC_MSG_CHECKING(for debugging environment)
AC_ARG_ENABLE(debug,[--enable-debug turn on debug compiler. default = no ] , [ AC_MSG_RESULT('Enabling debug') ] , [ AC_MSG_RESULT('Disabling debug')  ])
if test "x$enable_debug" = "xyes" ; then
	CXXFLAGS="-O0 -Wall -W -g3 -ggdb  -Werror "
	CFLAGS="-O0 -Wall -W -g3 -ggdb -p -Werror "

	AC_DEFINE(DEBUG,,[define if the DEBUG option is set])
else
dnl	AC_MSG_ERROR('Your debug mode is set to false')
	AC_MSG_RESULT('Your debug mode is set to false')
	AC_DEFINE(NDEBUG,,[no debug: assert not working])
fi
dnl AC_MSG_RESULT([$enable_debug])
AC_SUBST(CXXFLAGS)
AC_SUBST(CFLAGS)
])
