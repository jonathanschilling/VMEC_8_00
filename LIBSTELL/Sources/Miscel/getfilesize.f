      SUBROUTINE getfilesize (infile, filesize)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in)  :: infile
      INTEGER, INTENT(out)       :: filesize
      INTEGER :: iunit=100, istat
      LOGICAL :: lexist, lopen, lvalid
!DEC$ IF DEFINED (RISC)
     	INTEGER(4), EXTERNAL :: ftell_ 
!DEC$ ELSE
     	INTEGER(4), EXTERNAL :: ftell 
!DEC$ ENDIF
      
      filesize = -1
      INQUIRE(file=infile, exist=lexist)
      IF (.not.lexist) RETURN

      lvalid = .false.
      DO WHILE (.not.lvalid)
         INQUIRE(iunit, exist=lexist, opened=lopen)
         lvalid = lexist .and. (.not.lopen)
         IF (lvalid) EXIT
         iunit = iunit+1
      END DO

	OPEN (unit=iunit, file=infile, position='append',
     1      status='old', iostat=istat)
      IF (istat .ne. 0) RETURN
	
!DEC$ IF DEFINED (RISC)
	filesize = ftell_(iunit)
!DEC$ ELSE
	filesize = ftell(iunit)
!DEC$ ENDIF

      CLOSE (iunit)

      END SUBROUTINE getfilesize
      
