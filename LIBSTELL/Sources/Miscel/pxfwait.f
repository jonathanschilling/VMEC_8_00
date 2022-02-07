      SUBROUTINE pxfwait_g(istat, iretpid, ierror)
      IMPLICIT NONE
      INTEGER :: istat, iretpid, ierror
c !DEC$ IF DEFINED (WIN32)
c       istat = 0
c       iretpid = 0
c       ierror = 0
c !DEC$ ELSEIF .NOT.DEFINED (CRAY) .AND. .NOT. DEFINED(IRIX64)
      INTEGER, EXTERNAL :: wait

      iretpid = 0
      istat = 0

      ierror = wait(0)
c !DEC$ ELSE
c       CALL pxfwait(istat, iretpid, ierror)
c !DEC$ ENDIF

      END SUBROUTINE pxfwait_g
