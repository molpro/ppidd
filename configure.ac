dnl regnerate configure with: autoconf -Wall && rm -rf autom4te.cache
AC_INIT(PPIDD,master)
AC_PREREQ(2.69)
AC_LANG([C++])
AC_CONFIG_FILES([config.mk:src/configure/config.mk.in] [src/ppidd.h] [src/ppidd_module.F90])
AC_CONFIG_FILES([ppidd-config],[chmod +x ppidd-config])
AC_CONFIG_HEADERS([src/ppidd_config.h:src/configure/ppidd_config.h.in])
AC_CONFIG_AUX_DIR([src/configure])

AC_CONFIG_LINKS([GNUmakefile:GNUmakefile src/GNUmakefile:src/GNUmakefile])

AC_CANONICAL_BUILD()
AC_CANONICAL_HOST()

AC_SUBST([ppidd_srcdir],["`cd ${srcdir}; /bin/pwd`"])

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

AC_ARG_ENABLE(fortran,[AS_HELP_STRING([--enable-fortran],[build Fortran interface])])
AC_SUBST(enable_fortran)
AC_ARG_ENABLE(fortran-int64,[AS_HELP_STRING([--enable-fortran-int64],[specify Fortran integer size as 64-bit])],[fortran_integer="8"],[fortran_integer="4"])
AC_ARG_WITH(ga,[AS_HELP_STRING([--with-ga[=PATH]],[build using Global Arrays, PATH is the location of the GA include files])],,[with_ga="no"])
AC_ARG_ENABLE(nxtval,[AS_HELP_STRING([--disable-nxtval],[disable the helper process (PPIDD_Nxtval)])],,[enable_nxtval="yes"])
AC_ARG_WITH(mpi2,[AS_HELP_STRING([--with-mpi2],[build using MPI])],,[with_mpi2="no"])

AC_CHECK_TOOL([AR], [ar])
AC_SUBST(ARFLAGS,"-rS")
AC_PROG_RANLIB()
test "x${RANLIB}" = "x:" && RANLIB="true"
AC_PROG_INSTALL()

AC_PROG_CXX()
PPIDD_CXX_PP_TEST([if CXX is at least C++11 standard],[__cplusplus >= 201103L],,[AC_MSG_ERROR([require at least C++11 standard compiler])])

if test "x${enable_nxtval}" = xno ; then
 AC_MSG_NOTICE([Helper process will not be compiled; PPIDD_Nxtval will be disabled])
 AC_DEFINE(NO_NXTVAL_SERVER)
fi

if test "x${enable_fortran}" = xyes ; then
AC_PROG_FC()
m4_include([src/configure/may_fc_module_filename_format.m4])
MAY_FC_MODULE_FILENAME_FORMAT()
AC_FC_MODULE_FLAG()
else
FC=""
fi
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])

if test "x${ac_cv_sizeof_int}" = "x${fortran_integer}" ; then
 AC_DEFINE([FORTINT],[int])
elif test "x${ac_cv_sizeof_long}" = "x${fortran_integer}" ; then
 AC_DEFINE([FORTINT],[long])
elif test "x${ac_cv_sizeof_long_long}" = "x${fortran_integer}" ; then
 AC_DEFINE([FORTINT],[long long])
else
 AC_MSG_ERROR([Problem determining FORTINT])
fi

if test "x${with_ga}" != xno && test "x${with_mpi2}" != xno ; then
 AC_MSG_ERROR([both --with-ga and --with-mpi2 specified])
elif test "x${with_ga}" != xno ; then
 AC_MSG_NOTICE([PPIDD with Global Arrays])
 AC_SUBST(PPIDD_IMPL_DEFAULT,PPIDD_IMPL_GA_MPI)
elif test "x${with_mpi2}" != xno ; then
 AC_MSG_NOTICE([PPIDD with MPI2])
 AC_SUBST(PPIDD_IMPL_DEFAULT,PPIDD_IMPL_MPI2)
else
 AC_MSG_NOTICE([PPIDD without MPI])
 AC_SUBST(PPIDD_IMPL_DEFAULT,[PPIDD_IMPL_NO_MPI])
fi

if test "x${with_ga}" != xno ; then
 if test "x${with_ga}" != xyes ; then
  CPPFLAGS="${CPPFLAGS} -I${with_ga}"
 fi
 AC_CHECK_HEADERS([ga.h],,AC_MSG_ERROR([unable to locate ga.h header file]))
fi
if test "x${with_ga}" != xno || test "x${with_mpi2}" != xno ; then
 PPIDD_CXX_PP_TEST([if CXX is intel],[defined(__INTEL_COMPILER)],[CPPFLAGS="${CPPFLAGS} -DMPICH_SKIP_MPICXX"])
 AC_CHECK_FUNC([MPI_Init],,AC_MSG_ERROR([unable to locate MPI library]))
 AC_CHECK_HEADERS([mpi.h],,AC_MSG_ERROR([unable to locate mpi.h header file]))
fi

AC_SUBST_FILE(PPIDD_DEFINES_H)
PPIDD_DEFINES_H="${ppidd_srcdir}/src/ppidd_defines.h"
AC_SUBST_FILE(PPIDD_PROTOTYPES_H)
PPIDD_PROTOTYPES_H="${ppidd_srcdir}/src/ppidd_prototypes.h"

AC_OUTPUT()
