AC_PREREQ(2.59)
AC_INIT([himpi], [1.0.0], [http://hcl.ucd.ie/project/himpi])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE

AC_PROG_CC([mpicc])
AM_PROG_CC_C_O
AC_PROG_CXX([mpic++ mpicxx])
AC_ARG_ENABLE([smpi], [AC_HELP_STRING([--enable-smpi], [Enable the compilation with simgrid and smpi])], [smpi="yes"], [smpi="no"])
if test "$smpi" = "yes"; then
  CC=smpicc
  CXX=smpicc
  AC_DEFINE([HIMPI_SMPI], [1], [Smpi enabled])
fi

AC_PROG_LIBTOOL
AC_PROG_INSTALL

# Used to install python scripts with the specified python version
AM_PATH_PYTHON([2.6])


AC_MSG_CHECKING([whether debug is enabled])
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable debug (default: disabled)])], [debug="yes"], [debug="no"])
if test "$debug" = "yes"; then
	AC_DEFINE([HIMPI_DEBUG], [1], [Debug enabled])
	CFLAGS="-g -O0 -Wall"
	CXXFLAGS="-g -O0 -Wall -Wno-deprecated"
fi
AC_MSG_RESULT([$debug])


AC_MSG_CHECKING([whether use group config file])
AC_ARG_WITH([group_conf], [AC_HELP_STRING([--with-group-conf=FILE], [Group config file])])
if test "$with_group_conf" != ""; then
    if test -f "$with_group_conf"
	then
		AC_DEFINE_UNQUOTED([HIMPI_GROUP_CONFIG], ["$with_group_conf"], [Group config enabled])
	else
		AC_MSG_ERROR([The specified group config file does not exist.])
	fi
fi


AC_MSG_CHECKING([The HIMPI collectives that should be optimized by default])
AC_ARG_WITH([himpi_opid], [AC_HELP_STRING([--with-himpi-opid=<id>], [HiMPI operation/collective id])])
AC_DEFINE_UNQUOTED([HIMPI_OPID], ["$with_himpi_opid"], [Target HiMPI operation selected])


AC_MSG_CHECKING([Enable config generation])
AC_ARG_ENABLE([config_generation], [AC_HELP_STRING([--enable-config-generation], [Enable config-generation (default: disabled)])], [config_generation="yes"], [config_generation="no"])
if test "$config-generation" = "yes"; then
	AC_DEFINE([HIMPI_GENERATE_CONFIG], [1], [Config generation enabled])
fi
AC_MSG_RESULT([$config_generation])


AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_CHECK_LIB([gomp], [main])

ARCH=`uname -m`
if [[ $ARCH == "x86_64" ]];
then
  SUFFIX="64"
else
  SUFFIX=""
fi

AC_CHECK_LIB([rt], [main])
AC_CHECK_LIB([m], [main], [], [AC_MSG_ERROR([m required but not found])])
AC_CHECK_LIB([dl], [main], [], [AC_MSG_ERROR([dl required but not found])])
AC_DEFINE([_GNU_SOURCE], [1], [GNU])

AC_ARG_WITH([gsl_dir], [AC_HELP_STRING([--with-gsl-dir=DIR], [GNU Scientific Library directory])])
if test "$with_gsl_dir" != ""; then
	CPPFLAGS="$CPPFLAGS -I$with_gsl_dir/include"
	LDFLAGS="$LDFLAGS -L$with_gsl_dir/lib"
fi
AC_CHECK_LIB([gslcblas], [main], [], [AC_MSG_ERROR([gslcblas required but not found])])
AC_CHECK_LIB([gsl], [main], [], [AC_MSG_ERROR([gsl required but not found])])

AC_ARG_WITH([nlopt_dir], [AC_HELP_STRING([--with-nlopt-dir=DIR], [NLopt directory])])
if test "$with_nlopt_dir" != ""; then
	CPPFLAGS="$CPPFLAGS -I$with_nlopt_dir/include"
	LDFLAGS="$LDFLAGS -L$with_nlopt_dir/lib"
fi
if test "$smpi" = "no"; then
  AC_LANG_PUSH(C++)
  AC_CHECK_LIB([nlopt_cxx], [main], [], [nlopt_cxx="no"])
  AC_LANG_POP(C++)
else
  nlopt_cxx="no"
fi
AM_CONDITIONAL([HAVE_NLOPT_CXX], [test "$nlopt_cxx" != "no"])

# search pthread_attr_setaffinity_np within libpthread
saveCFLAGS="$CFLAGS"
saveLIBS="$LIBS"
CFLAGS="$CFLAGS $PTHREAD_CPPFLAGS $PTHREAD_CFLAGS"
LIBS="$LIBS $PTHREAD_LDFLAGS $PTHREAD_LIBS"

AC_CHECK_FUNCS([pthread_attr_setaffinity_np],
  [AC_DEFINE([USE_SCHED_AFFINITY],1,[Use affinity interface])])
CFLAGS="$saveCFLAGS"
LIBS="$saveLIBS"

#BG/P configuration
AC_CHECK_HEADER(mpix.h, [AC_DEFINE([HAVE_MPIX_H], [1], [Use MPI extensions for BG/P])])


# BLAS configuration
AC_MSG_CHECKING([for BLAS])
AC_ARG_WITH([blas], [AC_HELP_STRING([--with-blas=NAME], [BLAS (Basic Linear Algebra Subprograms): gsl, atlas, goto, essl, mkl(default: gsl)])])
if test "$with_blas" = ""; then
	with_blas="gsl"
fi
AC_MSG_RESULT([$with_blas])
case "$with_blas" in
gsl)
	AC_DEFINE([WITH_BLAS_GSL], [1], [Use the GSL BLAS])
	;;
atlas)
	AC_ARG_WITH([atlas_dir], [AC_HELP_STRING([--with-atlas-dir=DIR], [Automatically Tuned Linear Algebra Software directory])])
	if test "$with_atlas_dir" != ""; then
		CPPFLAGS="$CPPFLAGS -I$with_atlas_dir"
		LDFLAGS="$LDFLAGS -L$with_atlas_dir"
	fi
	AC_CHECK_LIB([f77blas], [main], [], [AC_MSG_ERROR([f77blas required but not found])])
	;;
goto)
	AC_ARG_WITH([goto_dir], [AC_HELP_STRING([--with-goto-dir=DIR], [Goto BLAS directory])])
	if test "$with_goto_dir" != ""; then
		CPPFLAGS="$CPPFLAGS -I$with_goto_dir/include"
		LDFLAGS="$LDFLAGS -L$with_goto_dir/lib"
	fi
	AC_CHECK_LIB([gfortran], [main])
	AC_CHECK_LIB([goto2], [main], [], [AC_MSG_ERROR([goto2 required but not found])])
	;;
essl)
        AC_ARG_WITH([essl_dir], [AC_HELP_STRING([--with-essl-dir=DIR], [ESSL directory])])
        if test "$with_essl_dir" != ""; then
                CPPFLAGS="$CPPFLAGS -I$with_essl_dir/include"
                LDFLAGS="$LDFLAGS -L$with_essl_dir/lib -L/opt/ibmcmp/xlf/bg/11.1/lib"
        fi
        AC_CHECK_LIB([xlomp_ser], [main], [], [AC_MSG_ERROR([xlomp_ser required but not found])])
        AC_CHECK_LIB([xlf90_r], [main], [], [AC_MSG_ERROR([xlf90_r required but not found])])
        AC_CHECK_LIB([xlfmath], [main], [], [AC_MSG_ERROR([xlfmath required but not found])])
        AC_CHECK_LIB([esslbg], [main], [], [AC_MSG_ERROR([esslbg required but not found])])
        # LIBS="-lmpitrace $LIBS"
        ;;
mkl)
	AC_ARG_WITH([mkl_dir], [AC_HELP_STRING([--with-mkl-dir=DIR], [Math Kernel Library directory])])
	if test "$with_mkl_dir" != ""; then
		CPPFLAGS="$CPPFLAGS -I$with_mkl_dir/include"
		LDFLAGS="$LDFLAGS -L$with_mkl_dir/lib"
	fi
	AC_CHECK_LIB([mkl_core], [main], [], [AC_MSG_ERROR([mkl_core required but not found])], [-lmkl_sequential]) 
	AC_CHECK_LIB([mkl_sequential], [main], [], [AC_MSG_ERROR([mkl_sequential required but not found])]) 
	AC_CHECK_LIB([mkl_gnu_thread], [main], [], [AC_MSG_ERROR([mkl_gnu_thread required but not found])])
	AC_CHECK_LIB([mkl_intel_lp64], [main], [], [AC_MSG_ERROR([mkl_intel_lp64 required but not found])]) 
	;;
*)
	AC_MSG_ERROR([not supported])
	;;
esac


AC_CONFIG_FILES([
	Makefile
	tools/Makefile
	himpi/Makefile
	tests/Makefile
	tests/mpi/Makefile
        tests/linear_algebra/Makefile
])

AC_CONFIG_SUBDIRS([MPIBlib])

DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN(himpi)

#define([svn_revision], esyscmd([sh -c "svnversion -n"]))
dnl
#DX_ENV_APPEND(SVN_REVISION, [svn_revision])




AC_OUTPUT
