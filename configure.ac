AC_INIT([intronProspector], 1.1.0)

##
# rebuild with:
#   autoreconf
#   autoconf
# requires autoconf-archive
##

PACKAGE_URL='https://github.com/diekhans/intronProspector'
AC_SUBST([PACKAGE_NAME])
AC_SUBST([PACKAGE_URL])
AC_SUBST([PACKAGE_VERSION])
AC_LANG([C++])
AC_PATH_PROG([PKG_CONFIG], [pkg-config], [no])

# Could not get ax_with_htslib.m4 from autoconf archive to work, so just hand
# coded it.
AC_ARG_WITH([htslib],
  [AS_HELP_STRING([--with-htslib=DIR],
    [use htlib installed in DIR (not the source tree), otherwise pkgconfig is check])], [])

AC_ARG_ENABLE([static-link],
              [AS_HELP_STRING([--enable-static-link],
                              [Link the program statically])],
              [case "${enableval}" in
               yes) static_link=true ;;
               no)  static_link=false ;;
               *) AC_MSG_ERROR([bad value ${enableval} for --enable-static-link]) ;;
               esac],
              [static_link=false])
AC_PROG_CC()
AC_PROG_CXX()

# FIXME: some of this is wonky, including using HTSLIB_EXTRA and
# having to add -rpath.  Need to dig into this

AC_MSG_NOTICE([checking for location of htslib include and library])

HTSLIB_EXTRA="-lcurl -lssl -lcrypt -lbz2 -lz -llzma"

if test x"${with_htslib}" != "x" ; then
   # --with-htslib=DIR specified
   if test "${with_htslib}" = "yes" ; then
       AC_MSG_ERROR([It appears that --with-htslib was supplied without an argument, please use --with-htslib=DIR])
   fi
   CPPFLAGS="${CPPFLAGS} -I${with_htslib}/include"
   AC_MSG_NOTICE(checking for htslib headers using ${CPPFLAGS})
   AC_CHECK_HEADER([htslib/sam.h],
       [AC_MSG_NOTICE([htslib/sam.h found])],
       [AC_MSG_ERROR([htslib/sam.h not found])])

   ldrpath="-L${with_htslib}/lib -Wl,-rpath,${with_htslib}/lib -Wl,-rpath-link,${with_htslib}/lib"
   LIBS="${LIBS} ${ldrpath} -lhts ${HTSLIB_EXTRA}"
   AC_MSG_NOTICE(checking for htslib library using ${LIBS})
   AC_CHECK_LIB(hts, hts_version,
      [AC_MSG_NOTICE([-lhts found])],
      [AC_MSG_ERROR([-lhts not found])])
else
    # use pkgconfig
    PKG_CHECK_MODULES([htslib], [htslib >= 0.20.0],
        [AC_MSG_RESULT([htslib found])],
        [AC_MSG_ERROR([pkg-config did not find htslib])])
    CPPFLAGS="${CPPFLAGS} ${htslib_CFLAGS}"

    # FIXME: with conda get warnings like
    # ld: warning: libdeflate.so.0, needed by .../lib/libhts.so, not found (try using -rpath or -rpath-link)

    htslibdir=`pkg-config --libs-only-L htslib | sed -e 's/^-L//' -e 's/ *$//'`
    ldrpath="-Wl,-rpath,${htslibdir} -Wl,-rpath-link,${htslibdir}"
    LIBS="${LIBS} ${ldrpath} ${htslib_LIBS} ${HTSLIB_EXTRA}"
fi

if test x$static_link = xtrue; then
    LDFLAGS="$LDFLAGS -static"
fi

AC_SUBST([CXXFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])
   
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile src/version.hh man/Makefile])

AC_OUTPUT
