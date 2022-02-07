      SUBROUTINE evolve(time_step, ier_flag, liter_flag, lscreen,
     1                  r0dot)
      USE vmec_main
      USE vmec_params, ONLY: bad_init_jac_flag
      USE vsvd
      USE xstuff
      USE precon2d
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: time_step, r0dot
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(inout) :: liter_flag
      LOGICAL, INTENT(in)  :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fsq1_threshold = 1.E-61_dp, 
     1                          r0dot_threshold = 5.E-06_dp
      CHARACTER*(*), PARAMETER :: qnewt_string =
     1      ' Turning on quasi-newton step at iteration: '
      REAL(rprec) :: fsq1, dtau, b1, bprec, fac
      LOGICAL :: lfinal_mesh
!      LOGICAL, SAVE :: lqnewt
C-----------------------------------------------
!
!     COMPUTE MHD FORCES
!
      lfinal_mesh = (ns .eq. MAXVAL(ns_array)) .and. .not.lprec2d

      IF (iter2 .lt. 100) THEN
         lprec2d = .false.
!         lqnewt = .false.
      ELSE IF (lfinal_mesh .and. 
     1        (fsqr1+fsqz1+fsql1).lt.fsq1_threshold
     2                     .and. r0dot.lt.r0dot_threshold) THEN  
         CALL compute_blocks(xc,xcdot,gc)

         time_step = 0.50_dp  !0.35_dp   !.25_dp
         iter1 = iter2-1; fsq = fsqr1 + fsqz1 + fsql1
!        nstep = 1
         xcdot = 0
      END IF

      CALL funct3d (lscreen, ier_flag)

!
!     COMPUTE ABSOLUTE STOPPING CRITERION
!
      IF (fsqr .le. ftolv .and. fsqz .le. ftolv .and.
     1    fsql .le. ftolv) THEN
         liter_flag = .false.
      ENDIF

      IF (irst.eq.2 .and. iter2.eq.1) ier_flag = bad_init_jac_flag

      IF (ier_flag.ne.0 .or. .not.liter_flag) RETURN

!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE

      IF (iter2 .eq. iter1) otau(:ndamp) = cp15/time_step

      IF (lprec2d) THEN
         bprec = 3
      ELSE
         bprec = 1
      END IF

      fsq1 = fsqr1 + fsqz1 + fsql1

      IF (iter2.gt.iter1 .and. fsq1.ne.zero) 
     1    dtau = MIN(ABS(1 - fsq1/fsq), bprec*cp15)   !MIN(ABS(LOG(fsq/fsq1)),cp15)

      fsq = fsq1

      IF (iter2 .le. 1) RETURN

      otau(1:ndamp-1) = otau(2:ndamp)

      IF (iter2 .gt. iter1) otau(ndamp) = dtau/time_step
      otav = SUM(otau(:ndamp))/ndamp
      dtau = time_step*otav/2

      b1  = one - dtau
      fac = one/(one + dtau)

!      IF (lprec2d .and. .not.lqnewt) THEN
!        IF ((fsqr + fsqz + fsql) .lt. fsq1_threshold) THEN
!           lqnewt = .true.
!           time_step = 1
!           DO iter = 1,2
!              IF (iter .eq. 1) iunit = 6
!              IF (iter .eq. 2) iunit = nthreed
!              WRITE (iunit,*) qnewt_string,iter2
!            END DO
!         END IF
!      END IF
!      IF (lqnewt) b1 = 0

!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN
!
      xcdot = fac*(xcdot*b1 + time_step*gc)
      xc    = xc + time_step*xcdot

      END SUBROUTINE evolve
