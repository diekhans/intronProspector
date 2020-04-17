AC_INIT([intronProspector], 1.0)
AC_LANG([C++])

# Could not get ax_with_htslib.m4 from autoconf archive to work, so just hand
# coded it.
AC_ARG_WITH([htslib],
  [AS_HELP_STRING([--with-htslib=DIR],
    [use the HTSlib source tree or installation in DIR. ])], [])

AC_PROG_CC()
AC_PROG_CXX()

AC_MSG_NOTICE([checking for location of HTSlib include and library])
hts_incl_found=no
hts_lib_found=no
if test x"${with_htslib}" != "x" ; then
   # explicit location
   hts_saved_CPPFLAGS=${CPPFLAGS}
   # try source tree layout
   HTSLIB_CPPFLAGS="-I${with_htslib}"
   CPPFLAGS="$CPPFLAGS $HTSLIB_CPPFLAGS"
   AC_MSG_NOTICE(checking for HTSlib header using $HTSLIB_CPPFLAGS)
   AC_CHECK_HEADER([htslib/sam.h],
       [hts_incl_found=yes],
       [HTSLIB_CPPFLAGS="-I${with_htslib}/include"
        CPPFLAGS="$CPPFLAGS $HTSLIB_CPPFLAGS"
        AC_MSG_NOTICE(checking for HTSlib header using $HTSLIB_CPPFLAGS)
        unset ac_cv_header_htslib_sam_h  # avoid cache
        AC_CHECK_HEADER([htslib/sam.h],
            [hts_incl_found=yes])])
   CPPFLAGS=${hts_saved_CPPFLAGS}
   # try source tree layout then lib dir
   hts_saved_LDFLAGS=${LDFLAGS}
   HTSLIB_LDFLAGS="-L${with_htslib} -Wl,-rpath,${with_htslib} -lhts"
   LDFLAGS="$LDFLAGS $HTSLIB_LDFLAGS"
   AC_MSG_NOTICE(checking for HTSlib library using $HTSLIB_LDFLAGS)
   AC_CHECK_LIB(hts, hts_version,
       [hts_lib_found=yes],
       [HTSLIB_LDFLAGS="-L${with_htslib}/lib -Wl,-rpath,${with_htslib} -lhts"
        LDFLAGS="$LDFLAGS $HTSLIB_LDFLAGS"
        unset ac_cv_lib_hts_hts_version  # avoid cache
        AC_MSG_NOTICE(checking for HTSlib library using $HTSLIB_LDFLAGS)
        AC_CHECK_LIB(hts, hts_version,
            [hts_lib_found=yes])])
   LDFLAGS=${hts_saved_LDFLAGS}
else
   # default location
   AC_MSG_NOTICE(checking for HTSlib header using standard include directory)
   AC_CHECK_HEADER([htslib/sam.h],
        [hts_incl_found=yes])
   AC_MSG_NOTICE(checking for HTSlib library using standard lib directory)
   AC_CHECK_LIB(hts, hts_version,
       [hts_lib_found=yes])
fi


if test "${hts_incl_found}" != "yes"  ; then
    AC_MSG_ERROR([HTSlib includes not found, use --with-htslib=DIR to select which HTSlib to use])
fi
if test "${hts_lib_found}" != "yes" ; then
    AC_MSG_ERROR([HTSlib library not found, use --with-htslib=DIR to select which HTSlib to use])
fi


CXXFLAGS="-O3 -g"
CPPFLAGS="${HTSLIB_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${HTSLIB_LDFLAGS} ${LDFLAGS}"
LIBS="-lhts -lcurl -lbz2 -lz -llzma"
AC_SUBST([CXXFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])
   
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])

AC_OUTPUT
