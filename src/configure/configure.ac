AC_INIT(QMay,master)
AC_CONFIG_FILES([config.mk:src/configure/config.mk.in])
AC_CONFIG_HEADERS([src/config.h:src/configure/config.h.in])
AC_CONFIG_AUX_DIR([src/configure])

AC_ARG_ENABLE(fortran,[AS_HELP_STRING([--enable-fortran],[build Fortran interface])])
AC_ARG_WITH(ga,[AS_HELP_STRING([--with-ga],[build using Global Arrays])],,[with_ga="no"])
AC_ARG_ENABLE(nxtval,[AS_HELP_STRING([--disable-nxtval],[disable the helper process (PPIDD_Nxtval)])],,[enable_nxtval="yes"])
AC_ARG_WITH(mpi2,[AS_HELP_STRING([--with-mpi2],[build using MPI])],,[with_mpi2="no"])

AC_SUBST(AR,"ar")dnl
AC_SUBST(ARFLAGS,"-rS")dnl
AC_PATH_PROG(AUTOCONF,autoconf)
AC_PATH_PROG(AUTOMAKE,automake)
AC_SUBST(AUTOMAKE_LIBDIR)dnl
if test "x${AUTOMAKE:+set}" = xset ; then
 if ${AUTOMAKE} --print-libdir 1>/dev/null 2>&1; then
  AUTOMAKE_LIBDIR="`${AUTOMAKE} --print-libdir`"
 fi
fi
AC_PATH_PROG(DOXYGEN,doxygen)
AC_PROG_RANLIB()
if test "x${RANLIB}" = "x:" ; then RANLIB="true"; fi

AC_PROG_CC()

if test "x${enable_nxtval}" = xno ; then
 AC_MSG_NOTICE([Helper process will not be compiled; PPIDD_Nxtval will be disabled])
 AC_DEFINE(NO_NXTVAL_SERVER)
fi

if test "x${enable_fortran}" = xyes ; then
AC_PROG_FC()
AC_ARG_ENABLE(fortran-integer,[AS_HELP_STRING([--enable-fortran-integer=4|8],[specify Fortran integer size])],,[enable_fortran_integer="4"])
case "x${enable_fortran_integer}" in
 x4 ) true ;;
 x8 ) AC_DEFINE(_I8_) ;;
 x* ) AC_MSG_ERROR([illegal value (${enable_fortran_integer}) given via --enable-fortran-integer option]) ;;
esac
fi

if test "x${with_ga}" != xno && test "x${with_mpi2}" != xno ; then
 AC_MSG_ERROR([both --with-ga and --with-mpi2 specified])
elif test "x${with_ga}" != xno ; then
 AC_MSG_NOTICE([PPIDD with Global Arrays])
 AC_DEFINE(GA_MPI)
elif test "x${with_mpi2}" != xno ; then
 AC_MSG_NOTICE([PPIDD with MPI2])
 AC_DEFINE(MPI2)
else
 true
 # serial
fi

AC_OUTPUT()
rm -f config.status config.log
