AC_INIT([aaphoto], [0.41], [mail@log69.com])
AM_INIT_AUTOMAKE([aaphoto], [0.41])
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# AC_PREREQ([2.67])
AC_CONFIG_SRCDIR([aaphoto.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgomp':
AC_CHECK_LIB([gomp], [main])
# FIXME: Replace `main' with a function in `-ljasper':
AC_CHECK_LIB([jasper], [main])
# FIXME: Replace `main' with a function in `-ljpeg':
AC_CHECK_LIB([jpeg], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpng':
AC_CHECK_LIB([png], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
