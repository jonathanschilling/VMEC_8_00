      SUBROUTINE getcarg(narg, arg, numargs)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)  :: narg
      INTEGER, INTENT(out) :: numargs
      CHARACTER*(*), INTENT(out) :: arg
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numchars
c !DEC$ IF DEFINED (LONESTAR) .OR. DEFINED (MCURIE)
c       INTEGER :: ier
c !DEC$ ENDIF
C-----------------------------------------------
c !DEC$ IF DEFINED (WIN32)
c       INTEGER nargs
c       numargs = nargs() - 1
c       CALL getarg(narg, arg, numchars)
c !DEC$ ELSEIF DEFINED (LINUX)
      INTEGER iargc
      numargs = iargc()
      CALL getarg(narg, arg)
c !DEC$ ELSEIF DEFINED (VMS)
c       CALL lib$get_foreign(arg,,numchars)
c       numargs = MIN(1,numchars)
c !DEC$ ELSEIF DEFINED (LONESTAR) .OR. DEFINED (MCURIE)
c       INTEGER iargc
c       numargs = iargc()
c       CALL pxfgetarg(narg, arg, numchars, ier)
c !DEC$ ELSE
c       INTEGER iargc, getarg
c       numargs = iargc()
c       numchars = getarg(narg, arg)
c !DEC$ ENDIF

      END SUBROUTINE getcarg
