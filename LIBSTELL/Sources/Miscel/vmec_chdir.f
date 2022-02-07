      INTEGER FUNCTION vmec_chdir(new_path)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in) :: new_path
c !DEC$ IF DEFINED (CRAY)
c       INTEGER :: ilen
c       iLEN = 0
c       CALL pxfchdir(new_path, ilen, vmec_chdir)
c !DEC$ ELSE
      INTEGER, EXTERNAL :: chdir

c !DEC$ IF DEFINED (SUNOS)
c       vmec_chdir = chdir(TRIM(new_path))
c !DEC$ ELSE
      vmec_chdir = chdir(TRIM(new_path) // char(0))
c !DEC$ ENDIF
c !DEC$ ENDIF
      END FUNCTION vmec_chdir
