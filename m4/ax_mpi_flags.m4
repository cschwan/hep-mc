AC_DEFUN([AX_MPI_FLAGS_INIT], [
    AC_LANG_CASE([C++],
                 [AC_ARG_VAR(MPICXX, [MPI C++ compiler command])
                  AC_PATH_PROGS(MPICXX, [mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC])])
])

AC_DEFUN([AX_MPI_FLAGS], [
    AC_LANG_CASE([C++], [
        _AX_MPI_FLAGS_GET_AND_CHECK(MPICXX, [AC_LANG_PROGRAM([[@%:@include <mpi.h>]],
            [[int r, s; MPI_Comm_rank(MPI_COMM_WORLD, &r); MPI_Comm_size(MPI_COMM_WORLD, &s);]])])
    ])
])

AC_DEFUN([_AX_MPI_FLAGS_GET_AND_CHECK], [
    AC_MSG_CHECKING([if $$1 supports -showme:compile])
    ax_mpi_flags=`$$1 -showme:compile 2>/dev/null`
    AS_IF([test "$?" -eq 0],
          [AC_MSG_RESULT([yes])
           ax_mpi_flags="${ax_mpi_flags} `$$1 -showme:link 2>/dev/null`"],
          [AC_MSG_RESULT([no])
           AC_MSG_CHECKING([if $$1 supports -showme])
           ax_mpi_flags=`$$1 -showme 2>/dev/null`
           AS_IF([test "$?" -eq 0],
                 [AC_MSG_RESULT([yes])],
                 [AC_MSG_RESULT([no])
                  AC_MSG_CHECKING([if $$1 supports -compile-info])
                  ax_mpi_flags=`$$1 -compile-info 2>/dev/null`
                  AS_IF([test "$?" -eq 0],
                        [AC_MSG_RESULT([yes])
                         ax_mpi_flags="${ax_mpi_flags} `$$1 -link-info 2>/dev/null`"],
                        [AC_MSG_RESULT([no])
                         AC_MSG_CHECKING([if $$1 supports -show])
                         ax_mpi_flags=`$$1 -show 2>/dev/null`
                         AS_IF([test "$?" -eq 0],
                               [AC_MSG_RESULT([yes])],
                               [AC_MSG_RESULT([no])])
                        ])
                 ])
          ])

    set dummy $ax_mpi_flags; shift
    for i in $[@]; do
        AS_CASE([$i],
                [-D*|-I*], [MPI_CPPFLAGS="$MPI_CPPFLAGS $i"],
                [-L*],     [MPI_LDFLAGS="$MPI_LDFLAGS $i"],
                [-l*],     [MPI_LIBS="$MPI_LIBS $i"],
                [echo "$as_me:$LINENO: ignored flag '$i'" >&AS_MESSAGE_LOG_FD])
    done

    MPI_CPPFLAGS=${MPI_CPPFLAGS# }
    MPI_LDFLAGS=${MPI_LDFLAGS# }
    MPI_LIBS=${MPI_LIBS# }

    OLD_CPPFLAGS=$CPPFLAGS
    OLD_LIBS=$LIBS
    OLD_LDFLAGS=$LDFLAGS

    CPPFLAGS=$MPI_CPPFLAGS
    LIBS=$MPI_LIBS
    LDFLAGS=$MPI_LDFLAGS

    AC_MSG_CHECKING([if extracted MPI flags work])
    AC_LINK_IFELSE([$2], [AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no])
        AC_MSG_FAILURE([failed to compile MPI test program.])])

    LDFLAGS=$OLD_LDFLAGS
    LIBS=$OLD_LIBS
    CPPFLAGS=$OLD_CPPFLAGS

    AC_SUBST(MPI_CPPFLAGS)
    AC_SUBST(MPI_LDFLAGS)
    AC_SUBST(MPI_LIBS)
])
