      SUBROUTINE vmec_getenv(ename, evalue)
      IMPLICIT NONE
      CHARACTER*(*) :: ename, evalue
c !DEC$ IF DEFINED (LONESTAR) .OR. DEFINED(MCURIE)
c       INTEGER :: lenname=0, lenval, ierror
c       CALL pxfgetenv(ename, lenname, evalue, lenval, ierror)
c !DEC$ ELSE
      CALL getenv(ename, evalue)
c !DEC$ ENDIF
      END SUBROUTINE vmec_getenv
