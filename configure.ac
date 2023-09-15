AC_INIT([KQED], [0.14], [renwick.james.hudspith@gmail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AM_PROG_AR
AC_C_BIGENDIAN
AC_CANONICAL_HOST

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

## Set up CUDA
AX_CUDA

## prefix!!
AC_DEFINE_UNQUOTED([HAVE_PREFIX],["$prefix"],[Look for where we prefix to])

## tell the code where the PRECOMP library is -> Here & only here!
AC_MSG_NOTICE([Setting LIB_PATH "$prefix"])
AC_DEFINE_UNQUOTED([LIB_PATH],["$prefix"],[Where the library and precomps are])
	  
## Tells us where we are prefixed to
AC_MSG_NOTICE([Configuring for location "$prefix"])
## If we prefix away from the default /usr/ then create binary GLU
## will still try and create shared object libraries there
AM_CONDITIONAL([PREF],[ test "$prefix" == "NONE" ])
   dnl if test "$prefix" == "NONE" ; then
   dnl     AC_MSG_NOTICE([Not creating binary "URFIT"])
   dnl else
   dnl     AC_MSG_NOTICE([Creating binary "URFIT" in $prefix/bin])
   dnl     mkdir -p $prefix/bin/
   dnl     cp -vp PRECOMP/REG* $prefix/bin
   dnl fi
AC_SUBST(PREF)

## C standard libraries
AC_CHECK_HEADERS([math.h],[LDFLAGS="$LDFLAGS"],
        [AC_MSG_ERROR([We need math.h])])

## C standard libraries
AC_CHECK_HEADERS([omp.h],[],[])

## Check for sys/time.h
AC_CHECK_HEADERS([sys/time.h],
	[AC_MSG_NOTICE([We are intrinsically timing our operations.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([time.h],
	[AC_MSG_NOTICE([We are providing an accurate date.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([immintrin.h],
	[AC_MSG_NOTICE([We are vectorising the crc32c checksums.])],[])

## I doubt we will not have these as they are all c-standard, check anyway
AC_CHECK_HEADERS([stdbool.h stdint.h strings.h stdio.h stdlib.h],[],[
AC_MSG_ERROR([Sorry, but we need these headers (strings.h , stdlib.h , stdio.h) somewhere in the code])])

## Intel MKL
AC_ARG_ENABLE([MKL],
	[ --enable-MKL],
	[case "${enableval}" in
	yes) AC_DEFINE([MKL],,[Link to the MKL library])
	     AC_MSG_NOTICE([MKL set])
	    LDFLAGS="$LDFLAGS -lmkl";;
	no) AC_MSG_NOTICE([MKL not set]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-MKL]) ;;
	esac],[])

## VERBOSE
AC_ARG_ENABLE([VERBOSE],
	[ --enable-VERBOSE],
	[case "${enableval}" in
	yes) AC_DEFINE([VERBOSE],,[Print out some more information than typical])
	     AC_MSG_NOTICE([VERBOSE set]) ;;	
	no) AC_MSG_NOTICE([VERBOSE not set]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-VERBOSE]) ;;
	esac],[])

## Can't get enough of that wonderful DUFF
AC_ARG_ENABLE([DUFF],
	[ --enable-DUFF],
	[case "${enableval}" in
	yes) AC_DEFINE([DUFF],,[Duff's device loop unrolling])
	     AC_MSG_NOTICE([DUFF set]) ;;	
	no) AC_MSG_NOTICE([DUFF not set]) ;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-DUFF]) ;;
	esac],[])

## Some compiler checks
## My code is littered with consts
AC_C_CONST([])

## Set up config.h, one of the most important headers of them all
AC_CONFIG_HEADERS([config.h])

## where are my files?
AC_CONFIG_FILES([
Makefile
src/Makefile
src_cuda/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
