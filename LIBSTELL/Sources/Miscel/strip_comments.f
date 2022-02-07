      SUBROUTINE strip_comments(input_file)
      USE safe_open_mod
      IMPLICIT NONE
!
!     strips comment lines (starting with '!') from input_file
!     renames clean file input_file // '.stripped'
!
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER*(*) :: input_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_strip = 20
      INTEGER :: istat, iustrip, iunew
      CHARACTER*(200) :: line
      LOGICAL :: lex
C-----------------------------------------------
      INQUIRE (file=input_file, exist = lex)
      IF (.not. lex) THEN
         PRINT *,TRIM(input_file),' does not exist!'
         STOP
      END IF
      iustrip = unit_strip
      CALL safe_open(iustrip, istat, input_file, 'old','formatted')
      IF (istat .ne. 0) THEN
        line = 'error opening ' // TRIM(input_file) //
     1         ' in STRIP_COMMENTS'
        PRINT *, line
        PRINT *,'istat = ', istat
        STOP
      END IF
      iunew = iustrip + 1
      CALL safe_open(iunew, istat, TRIM(input_file) // '.stripped',
     1     'replace', 'formatted')
      IF (istat .ne. 0) THEN
        line = 'error opening ' // TRIM(input_file) //
     1      '.stripped  in STRIP_COMMENTS'
        PRINT *, line
        PRINT *,'istat = ', istat
        STOP
      END IF
      DO
         READ(iustrip, '(a)', END=100) line
         line = ADJUSTL(line)
         IF (line(1:1) == '!') CYCLE
         WRITE(iunew, '(a)') TRIM(line)
      END DO

 100  CONTINUE

      CLOSE(iustrip)
      CLOSE(iunew)

      END SUBROUTINE strip_comments
