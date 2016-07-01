# ================================================================
#
# ================================================================
#
# SYNOPSIS
#
#   AX_MPI_FLAGS(MPI-COMPILER)
#
# DESCRIPTION
#
#    Uses MPI-COMPILER to query which flags need to be set to properly
#    compile MPI programs. This macro calls
#
#      AC_SUBST(MPI_CPPFLAGS)
#      AC_SUBST(MPI_LDFLAGS)
#      AC_SUBST(MPI_LIBS)
#
#    which are the preprocessor- and linker-flags and libraries needed for
#    MPI programs. If MPI-COMPILER does not exist or this macro is unable to
#    query the flags, it aborts configure with an error message.
#
#    In contrast to AX_MPI this macro allows to easily mix MPI- and non-MPI
#    programs.
#
# LICENSE
#
#   Copyright (c) 2014-2016 Christopher Schwan <handgranaten-herbert@posteo.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([_AX_MPI_FLAGS], [dnl
    AC_MSG_CHECKING([how to query MPI flags from $$1])

    AS_IF([ax_mpi_flags=`$$1 -showme:compile 2>/dev/null`; test "$?" -eq 0],
          [ax_mpi_flags="$ax_mpi_flags `$$1 -showme:link 2>/dev/null`"
           AC_MSG_RESULT([-showme:compile -showme:link])],

          [ax_mpi_flags=`$$1 -showme 2>/dev/null`; test "$?" -eq 0],
          [AC_MSG_RESULT([-showme])],

          [ax_mpi_flags=`$$1 -compile-info 2>/dev/null`; test "$?" -eq 0],
          [ax_mpi_flags="$ax_mpi_flags `$$1 -link-info 2>/dev/null`"
           AC_MSG_RESULT([-compile-info -link-info])],

          [ax_mpi_flags=`$$1 -show 2>/dev/null`; test "$?" -eq 0],
          [AC_MSG_RESULT([-show])],

          [AC_MSG_RESULT([])
           AC_MSG_ERROR([cannot extract MPI flags from $$1])]
    )

    AS_ECHO(["$as_me:$LINENO: extracted flags are $ax_mpi_flags"]) >&AS_MESSAGE_LOG_FD

    set dummy ${ax_mpi_flags}; shift
    while test $[@%:@] -gt 0; do
        AS_CASE([$[1]],
                [-D|-I],    [MPI_CPPFLAGS="$MPI_CPPFLAGS $[1] $[2]"; shift],
                [-D*|-I*],  [MPI_CPPFLAGS="$MPI_CPPFLAGS $[1]"],
                [-L],       [MPI_LDFLAGS="$MPI_LDFLAGS $[1] $[2]"; shift],
                [-L*|-Wl*], [MPI_LDFLAGS="$MPI_LDFLAGS $[1]"],
                [-l],       [MPI_LIBS="$MPI_LIBS $[1] $[2]"; shift],
                [-l*],      [MPI_LIBS="$MPI_LIBS $[1]"],

                [AS_ECHO(["$as_me:$LINENO: ignored flag '$[1]'"]) >&AS_MESSAGE_LOG_FD]
        )
        shift
    done

    MPI_CPPFLAGS=${MPI_CPPFLAGS## }
    MPI_LDFLAGS=${MPI_LDFLAGS## }
    MPI_LIBS=${MPI_LIBS## }

    OLD_CPPFLAGS=$CPPFLAGS
    OLD_LIBS=$LIBS
    OLD_LDFLAGS=$LDFLAGS

    CPPFLAGS=$MPI_CPPFLAGS
    LIBS=$MPI_LIBS
    LDFLAGS=$MPI_LDFLAGS

    AC_MSG_CHECKING([if extracted MPI flags work])
    AC_LINK_IFELSE([$2], [AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no])
        AC_MSG_FAILURE([cannot compile MPI test program])])

    LDFLAGS=$OLD_LDFLAGS
    LIBS=$OLD_LIBS
    CPPFLAGS=$OLD_CPPFLAGS

    AC_SUBST(MPI_CPPFLAGS)
    AC_SUBST(MPI_LDFLAGS)
    AC_SUBST(MPI_LIBS)
])

m4_define([_AX_MPI_FLAGS_C_TESTBODY],
          [AC_LANG_PROGRAM([[@%:@include <mpi.h>]],
           [[int r, s; MPI_Comm_rank(MPI_COMM_WORLD, &r); MPI_Comm_size(MPI_COMM_WORLD, &s)]])])

AC_DEFUN([AX_MPI_FLAGS], [
    AC_LANG_CASE([C++], [_AX_MPI_FLAGS([$1], [_AX_MPI_FLAGS_C_TESTBODY])])
])
