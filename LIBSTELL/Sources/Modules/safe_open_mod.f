      MODULE safe_open_mod
!
!     Module for performing a "safe" open of a file for
!     a Fortran read/write operation. Makes sure the requested file
!     unit number is not in use, and increments it until an unused
!     unit is found
!
      CONTAINS

      SUBROUTINE safe_open(iunit, istat, filename, filestat,
     1           fileform, record_in, access_in)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: iunit
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in), OPTIONAL :: record_in
      CHARACTER*(*), INTENT(in) :: filename, filestat, fileform
      CHARACTER*(*), INTENT(in), OPTIONAL :: access_in
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      CHARACTER*(*), PARAMETER :: cdelim = "apostrophe", 
     1     cform="formatted", cunform="unformatted", 
     2     cscratch="scratch", cseq="sequential"
      CHARACTER*(10) :: acc_type
      LOGICAL :: lopen, lexist, linvalid
C-----------------------------------------------
!
!     Check that unit is not already opened
!
      linvalid = .true.
      DO WHILE (linvalid)
         INQUIRE(iunit, exist=lexist, opened=lopen)
         linvalid = (.not.lexist) .or. lopen
         IF (.not.linvalid) EXIT
         iunit = iunit+1
      END DO

      IF (PRESENT(access_in)) THEN
         acc_type = TRIM(access_in)
      ELSE
         acc_type = cseq
      END IF

      lexist = (filestat(1:1).eq.'s') .or. (filestat(1:1).eq.'S')        !Scratch file

      IF (PRESENT(access_in)) THEN
         acc_type = TRIM(access_in)
      ELSE
         acc_type = 'SEQUENTIAL'
      END IF

      SELECT CASE (fileform(1:1))
      CASE ('u', 'U')
         IF (PRESENT(record_in)) THEN
            IF (lexist) THEN
               OPEN(unit=iunit, form=cunform, status=cscratch,
     1              recl=record_in, access=acc_type, iostat=istat)
            ELSE
               OPEN(unit=iunit, file=TRIM(filename), form=cunform,
     1              status=TRIM(filestat), recl=record_in,
     2              access=acc_type, iostat=istat)
            END IF
         ELSE
            IF (lexist) THEN
               OPEN(unit=iunit, form=cunform, status=cscratch,
     1              access=acc_type, iostat=istat)
            ELSE
               OPEN(unit=iunit, file=TRIM(filename), form=cunform,
     1              status=TRIM(filestat), access=acc_type,iostat=istat)
            END IF
         END IF

      CASE DEFAULT
         IF (PRESENT(record_in)) THEN
            IF (lexist) THEN
               OPEN(unit=iunit, form=cform, status=cscratch,
     1              delim=TRIM(cdelim), recl=record_in, access=acc_type,
     2              iostat=istat)
            ELSE
               OPEN(unit=iunit, file=TRIM(filename), form=cform,
     1              status=TRIM(filestat), delim=TRIM(cdelim),
     2              recl=record_in, access=acc_type, iostat=istat)
            END IF
         ELSE
            IF (lexist) THEN
               OPEN(unit=iunit, form=cform, status=cscratch,
     1              delim=TRIM(cdelim), access=acc_type, iostat=istat)
            ELSE
               OPEN(unit=iunit, file=TRIM(filename), form=cform,
     1             status=TRIM(filestat), delim=TRIM(cdelim),
     2             access=acc_type, iostat=istat)
            END IF
         END IF

      END SELECT

      END SUBROUTINE safe_open

      END MODULE safe_open_mod
