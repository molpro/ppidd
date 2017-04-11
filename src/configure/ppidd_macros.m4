dnl msg,test,[true],[false]
AC_DEFUN([PPIDD_CXX_PP_TEST],[
AC_MSG_CHECKING([$1])
AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[#if $2
#else
      choke me
#endif]])],
[AC_MSG_RESULT([yes])
$3],
[AC_MSG_RESULT([no])
$4])
AC_LANG_POP([C++])
])

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
if test "x${FC_VENDOR:+set}" != xset ; then
 AC_MSG_WARN([unable to determine FC_VENDOR])
fi
])
