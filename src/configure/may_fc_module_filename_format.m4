# SYNOPSIS
#
#   MAY_FC_MODULE_FILENAME_FORMAT
#
# DESCRIPTION
#
#   Find the Fortran 90 module filename format. The format is stored in
#   the variable FC_MODFMT and empty if it cannot be determined. The
#   following modifiers may appear in FC_MODFMT:
#
#      %m - lower case module name
#      %M - upper case module name
#
#   The result or "unknown" is cached in the cache variable
#   may_cv_fc_module_fmt.
#
# LICENSE
#
#   Copyright (c) 2014 Andy May <ajmay81@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([MAY_FC_MODULE_FILENAME_FORMAT],
[AC_CACHE_CHECK([Fortran 90 module filename format], [may_cv_fc_module_fmt],
[AC_REQUIRE([AC_FC_MODULE_EXTENSION]) dnl to set FC_MODEXT
AC_LANG_PUSH(Fortran)
mkdir conftest.dir
cd conftest.dir
may_cv_fc_module_fmt="unknown"
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'hello world'
      end subroutine conftest_routine
      end module conftest_module]],
  [if test -r "conftest_module.${FC_MODEXT}" ; then
    may_cv_fc_module_fmt="%m.${FC_MODEXT}"
   elif test -r "CONFTEST_MODULE.${FC_MODEXT}" ; then
    may_cv_fc_module_fmt="%M.${FC_MODEXT}"
   fi])
cd ..
rm -rf conftest.dir
AC_LANG_POP(Fortran)
])
FC_MODFMT="${may_cv_fc_module_fmt}"
if test "x${FC_MODFMT}" = xunknown; then
  FC_MODFMT=""
fi
AC_SUBST([FC_MODFMT])dnl
])
