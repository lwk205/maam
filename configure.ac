# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#autoscan genera una versione di configure.ac osservando
#la struttura del programma
#AC_DEFINE passa variabili al compilatore

AC_INIT([maam],
        [1.0],
        [Simone Furini sfurini@deis.unibo.it],
        [maam])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
if test "x$GCC" = "xyes" ; then
	CFLAGS="-g -Wall -O3 -funroll-loops"
	CXXFLAGS="-g -Wall -O3 -funroll-loops"
fi

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strchr strrchr feenableexcept])

######################################################################
#   conditional compiling                                            #
######################################################################
AC_ARG_ENABLE(forcemomenta,
		[  --enable-forcemomenta    compile b_forcemomenta],
		[case "${enableval}" in
   			yes) forcemomenta=true ;;
   			no)  forcemomenta=false ;;
       			*) AC_MSG_ERROR(bad value ${enableval} for --enable-forcemomenta) ;;
	     	esac],[forcemomenta=false])
AM_CONDITIONAL(FORCEMOMENTA, test x$forcemomenta = xtrue)


AC_CONFIG_FILES([
   Makefile
   README
   m4/Makefile
   doc/Makefile
   dat/Makefile
   include/Makefile
   lib/Makefile
   test/Makefile
   src/Makefile
   src/laos/Makefile
   src/nepal/Makefile
])

AC_OUTPUT
