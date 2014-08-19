AC_INIT(QMay,master)
AC_CONFIG_FILES([config.mk:src/configure/config.mk.in])
AC_CONFIG_AUX_DIR([src/configure])

AC_PATH_PROG(AUTOCONF,autoconf)
AC_PATH_PROG(AUTOMAKE,automake)
AC_SUBST(AUTOMAKE_LIBDIR)dnl
if test "x${AUTOMAKE:+set}" = xset ; then
 if ${AUTOMAKE} --print-libdir 1>/dev/null 2>&1; then
  AUTOMAKE_LIBDIR="`${AUTOMAKE} --print-libdir`"
 fi
fi
AC_PATH_PROG(DOXYGEN,doxygen)

AC_OUTPUT()
rm -f config.status config.log
