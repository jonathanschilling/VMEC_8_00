      SUBROUTINE pxffork_g (ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
c !DEC$ IF DEFINED (WIN32)
c       ierror = 0
c       ipid = 0
c !DEC$ ELSEIF .NOT.DEFINED (CRAY) .AND. .NOT. DEFINED(IRIX64)
      INTEGER, EXTERNAL :: fork

      ierror = 0
      ipid = fork()
      IF (ipid < 0) ierror = -ipid
c !DEC$ ELSE
c       CALL pxffork (ipid, ierror)
c !DEC$ ENDIF
      END SUBROUTINE pxffork_g
