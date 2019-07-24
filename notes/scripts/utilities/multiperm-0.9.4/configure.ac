#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([multiperm], [0.9.4], [anandam@u.washington.edu])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([options.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lR':
AC_CHECK_LIB([R], [main])
# FIXME: Replace `main' with a function in `-lRblas':
AC_CHECK_LIB([Rblas], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_TYPE_SIGNAL
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor setenv sqrt strcasecmp strrchr])

# Enable or disable the use of R
AC_ARG_ENABLE(R,
	AC_HELP_STRING([--enable-R], [Use the R statistical package for random number generation and pvalue calculations]),
	[use_r=${enableval}],
	[use_r=no])
AM_CONDITIONAL([USE_R], [test "$use_r" = "yes"])

# Find out where R is (by default, it's /usr/local/lib/R)
# Make sure your PATH is set up to run R
AC_CHECK_PROG([RHOME], [R], [`R RHOME`])

# Enable or disable the use the rnazRandomizeAln.pl script for comparisons
AC_ARG_ENABLE(compRNAz,
	AC_HELP_STRING([--enable-compRNAz], [Enable the use of the rnazRandomizeAln.pl script for comparisons with multiperm]),
	[use_rnaz=${enableval}],
	[use_rnaz=no])
AM_CONDITIONAL([USE_RNAZ], [test "$use_rnaz" = "yes"])

# Find out where rnazRandomizeAln.pl is
# Make sure your PATH is set up to run rnazRandomizeAln.pl
AC_CHECK_PROG([RNAZ_RANDOMIZE_ALN], [rnazRandomizeAln.pl], [`which rnazRandomizeAln.pl`])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
