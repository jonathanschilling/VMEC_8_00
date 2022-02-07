      SUBROUTINE runvmec(input_file, iseq_count, lreseta, ier_flag,
     1   ireset, lfirst_call, lscreen, reset_file_name)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, more_iter_flag
      USE vsvd
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iseq_count, ier_flag, ireset
      LOGICAL :: lreseta, lfirst_call, lscreen
      CHARACTER*(*) :: input_file, reset_file_name
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: igrid0, igrid, nsval, index_end, index_dat,
     1   ns_min, ier_init, jacob_off
      LOGICAL :: interp_flag, lreset
C-----------------------------------------------

!
!
!     INDEX OF LOCAL VARIABLES
!
!     ier_flag   specifies error condition IF nonzero
!     interp_flag
!                = FALSE,compute xc array from scaled boundary data
!                = TRUE, compute xc array by interpolation
!     lfirst_call
!                = T, initial CALL to runvmec
!                = F, runvmec has been previously called
!
      CALL second0 (timeon)

!
!     PARSE input_file into path/input.ext
!
      index_dat = index(input_file,'input.')
      index_end = len_trim(input_file)
      IF (index_dat .gt. 0) THEN
         input_extension  = input_file(index_dat+6:index_end)
      ELSE
         input_extension = input_file(1:index_end)
         input_file = 'input.'//TRIM(input_extension)
      END IF

!
!     INITIALIZE PARAMETERS
!
      lreset = .false.
      ier_init = ier_flag
      IF (lfirst_call .or. lreseta) lreset = .true.

      IF (ier_init.ne.more_iter_flag .and. 
     1    ier_init.ne.bad_jacobian_flag) THEN
         CALL vsetup (lreset, iseq_count)
      ELSE
         iequi = 0
         IF (lfreeb) ivac = 1    !!Must restart vacuum calculations IF free boundary
      END IF

!
!     READ INPUT FILES INDATA, MGRID_FILE
!     USE ISEQ_COUNT TO AVOID REREADING MGRID FILE
!
      CALL readin
     1    (input_file, iseq_count, lfirst_call, ier_flag, lscreen)
      IF (ier_flag .ne. 0) GOTO 1000

!
!     COMPUTE INVARIANT ARRAYS
!
      CALL fixaray

      IF (ier_init .ne. 4) THEN
         WRITE (nthreed, 230)
  230 FORMAT(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)
      END IF

!
!     COMPUTE INITIAL SOLUTION ON COARSE GRID
!     IF PREVIOUS SEQUENCE DID NOT CONVERGE, DO COARSE RESTART
!

      IF (lreseta) THEN
        igrid0 = 1
      ELSE
        igrid0 = multi_ns_grid
      ENDIF


      imovephi = 0
      ns_min = 0
      jacob_off = 0
      ier_flag = ier_init
      IF (ier_flag .eq. bad_jacobian_flag) jacob_off = 1
      IF (ALL(ns_array .eq. 0)) THEN
         ier_flag = 8
         GOTO 1000
      END IF

      DO igrid = igrid0, multi_ns_grid + jacob_off
         IF (jacob_off.eq.1 .and. igrid.eq.igrid0) THEN
!           TRY TO GET NON-ZERO JACOBIAN ON A 3 PT RADIAL MESH
            nsval = 3
            ftolv = 1.e-7_dp
         ELSE
            nsval = ns_array(igrid-jacob_off)
            IF (nsval .le. ns_min) CYCLE
            ns_min = nsval
            ftolv = ftol_array(igrid-jacob_off)
         END IF
         IF (igrid .eq. igrid0) THEN
            interp_flag = .false.
         ELSE
            interp_flag = .true.
         ENDIF
         CALL eqsolve (nsval, interp_flag, ier_flag, lreseta,
     1                 lfirst_call, lscreen, reset_file_name)
         IF (imatch_phiedge .eq. 3) imovephi = 1
         IF (ier_flag .ne. 0 .and. ier_flag .ne. 4) EXIT
      END DO

  100 CONTINUE

      CALL second0 (timeoff)
      timer(tsum) = timer(tsum) + timeoff - timeon

!
!     WRITE OUTPUT TO THREED1, WOUT FILES; FREE MEMORY ALLOCATED GLOBALLY
!
 1000 CALL fileout (iseq_count, ier_flag, ireset, lscreen)

      END SUBROUTINE runvmec
