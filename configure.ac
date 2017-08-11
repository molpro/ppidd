dnl regnerate configure with: autoconf -Wall && rm -rf autom4te.cache
AC_INIT(PPIDD,master)
AC_PREREQ(2.69)
AC_LANG([C++])
AC_CONFIG_FILES([config.mk:src/configure/config.mk.in])
AC_CONFIG_HEADERS([src/ppidd_config.h:src/configure/ppidd_config.h.in])
AC_CONFIG_AUX_DIR([src/configure])

AC_CANONICAL_BUILD()
AC_CANONICAL_HOST()

dnl msg,test,[true],[false]
AC_DEFUN([PPIDD_CXX_PP_TEST],[
AC_MSG_CHECKING([$1])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[#if $2
#else
      choke me
#endif]])],
[AC_MSG_RESULT([yes])
$3],
[AC_MSG_RESULT([no])
$4])])

dnl msg,test,[true],[false]
AC_DEFUN([PPIDD_FC_PP_TEST],[
AC_MSG_CHECKING([$1])
AC_LANG_PUSH([Fortran])
AC_LANG_CASE([Fortran],[ac_ext=F])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[#if $2
#else
      choke me
#endif]])],
[AC_MSG_RESULT([yes])
$3],
[AC_MSG_RESULT([no])
$4])
AC_LANG_POP([Fortran])
])

AC_DEFUN([PPIDD_FC_VENDOR],[
for i in amd intel pgi pathscale sun g95 ibm gnu cray; do
 if test "x${FC_VENDOR:+set}" = xset ; then break; fi
 case "x${i}" in
  xamd       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__OPEN64__)],[FC_VENDOR="${i}"]) ;;
  xcray      ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(_CRAYFTN)],[FC_VENDOR="${i}"]) ;;
  xg95       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__G95__)],[FC_VENDOR="${i}"]) ;;
  xgnu       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__GFORTRAN__)],[FC_VENDOR="${i}"]) ;;
  xibm       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__xlc__) || defined(__xlC__)],[FC_VENDOR="${i}"]) ;;
  xintel     ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__INTEL_COMPILER)],[FC_VENDOR="${i}"]) ;;
  xpathscale ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__PATHSCALE__)],[FC_VENDOR="${i}"]) ;;
  xpgi       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__PGI)],[FC_VENDOR="${i}"]) ;;
  xsun       ) PPIDD_FC_PP_TEST([if FC is ${i}],[defined(__SUNPRO_F90)],[FC_VENDOR="${i}"]) ;;
  *          ) AC_MSG_ERROR([unknown vendor '${i}']) ;;
 esac
done
test "x${FC_VENDOR:+set}" = xset || AC_MSG_WARN([unable to determine FC_VENDOR])
])

AC_ARG_ENABLE(fortran,[AS_HELP_STRING([--enable-fortran],[build Fortran interface])])
AC_ARG_WITH(ga,[AS_HELP_STRING([--with-ga[=PATH]],[build using Global Arrays, PATH is the location of the GA include files])],,[with_ga="no"])
AC_ARG_ENABLE(nxtval,[AS_HELP_STRING([--disable-nxtval],[disable the helper process (PPIDD_Nxtval)])],,[enable_nxtval="yes"])
AC_ARG_WITH(mpi2,[AS_HELP_STRING([--with-mpi2],[build using MPI])],,[with_mpi2="no"])

AC_SUBST(AR,"ar")dnl
AC_SUBST(ARFLAGS,"-rS")dnl
AC_PATH_PROG(DOXYGEN,doxygen)
AC_PROG_RANLIB()
test "x${RANLIB}" = "x:" && RANLIB="true"
AC_PROG_INSTALL()

AC_PROG_CXX()
AC_PATH_PROG(CXX,${CXX})

if test "x${enable_nxtval}" = xno ; then
 AC_MSG_NOTICE([Helper process will not be compiled; PPIDD_Nxtval will be disabled])
 AC_DEFINE(NO_NXTVAL_SERVER)
fi

if test "x${enable_fortran}" = xyes ; then
AC_PROG_FC()
AC_PATH_PROG(FC,${FC})
m4_include([src/configure/may_fc_module_filename_format.m4])
MAY_FC_MODULE_FILENAME_FORMAT()
PPIDD_FC_VENDOR()
AC_ARG_ENABLE(fortran-integer,[AS_HELP_STRING([--enable-fortran-integer=4|8],[specify Fortran integer size])],,[enable_fortran_integer="4"])
case "x${enable_fortran_integer}" in
 x4 ) true ;;
 x8 ) case "x${FC_VENDOR}" in
       xamd       ) FCFLAGS="${FCFLAGS} -i8" ;;
       xcray      ) FCFLAGS="${FCFLAGS} -s integer64" ;;
       xg95       ) FCFLAGS="${FCFLAGS} -i8" ;;
       xgnu       ) FCFLAGS="${FCFLAGS} -fdefault-integer-8" ;;
       xibm       ) FCFLAGS="${FCFLAGS} -qintsize=8" ;;
       xintel     ) FCFLAGS="${FCFLAGS} -i8" ;;
       xpathscale ) FCFLAGS="${FCFLAGS} -i8" ;;
       xpgi       ) FCFLAGS="${FCFLAGS} -i8" ;;
       xsun       ) FCFLAGS="${FCFLAGS} -xtypemap=integer:64" ;;
       x*         ) AC_MSG_WARN([unable to determine flag for 8-byte Fortran integers, check FCFLAGS set correctly]);;
      esac ;;
 x* ) AC_MSG_ERROR([illegal value (${enable_fortran_integer}) given via --enable-fortran-integer option]) ;;
esac
else
dnl assume Fortran standard if compiling without Fortran compiler
enable_fortran_integer="4"
FC=""
fi
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])

if test "x${ac_cv_sizeof_int}" = "x${enable_fortran_integer}" ; then
 AC_DEFINE([FORTINT],[int])
 AC_DEFINE([FORTINT_MPI],[MPI_INT])
elif test "x${ac_cv_sizeof_long}" = "x${enable_fortran_integer}" ; then
 AC_DEFINE([FORTINT],[long])
 AC_DEFINE([FORTINT_MPI],[MPI_LONG])
elif test "x${ac_cv_sizeof_long_long}" = "x${enable_fortran_integer}" ; then
 AC_DEFINE([FORTINT],[long long])
 AC_DEFINE([FORTINT_MPI],[MPI_LONG_LONG])
else
 AC_MSG_ERROR([Problem determining FORTINT])
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
 AC_MSG_NOTICE([PPIDD for serial])
fi

if test "x${with_ga}" != xno ; then
 if test "x${with_ga}" != xyes ; then
  CPPFLAGS="${CPPFLAGS} -I${with_ga}"
 fi
 AC_CHECK_HEADER([ga.h],,AC_MSG_ERROR([unable to locate ga.h header file]))
fi
if test "x${with_ga}" != xno || test "x${with_mpi2}" != xno ; then
 PPIDD_CXX_PP_TEST([if CXX is intel],[defined(__INTEL_COMPILER)],[CPPFLAGS="${CPPFLAGS} -DMPICH_SKIP_MPICXX"])
 AC_CHECK_FUNC([MPI_Init],,AC_MSG_ERROR([unable to locate MPI library]))
 AC_CHECK_HEADER([mpi.h],,AC_MSG_ERROR([unable to locate mpi.h header file]))
fi

AC_OUTPUT()
