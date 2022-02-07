      SUBROUTINE vmec_system(cmd, ierror)
      INTEGER, OPTIONAL :: ierror
      INTEGER :: ireturn
      CHARACTER*(*), INTENT(in) :: cmd

c !DEC$ IF DEFINED (CRAY)
c       INTEGER, EXTERNAL :: ishell
c       ireturn = ishell(cmd)
c !DEC$ ELSEIF DEFINED (RISC)
c       CALL system(cmd, ireturn)
c !DEC$ ELSEIF DEFINED (LINUX) .OR. DEFINED(OSF1)
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd))
c !DEC$ ELSEIF DEFINED(WIN32) .OR. DEFINED(SUNOS)
c       INTEGER, EXTERNAL :: system
c       ireturn = system(TRIM(cmd))
c !DEC$ ELSE
c       INTEGER, EXTERNAL :: system
c       ireturn = system(TRIM(cmd) // char(0))
c !DEC$ ENDIF
      IF (PRESENT(ierror)) ierror = ireturn

      END SUBROUTINE vmec_system
