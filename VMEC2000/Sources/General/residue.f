      SUBROUTINE residue(gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rsc, rcs, rcc, meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE vsvd
      USE xstuff
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax) ::
     1  gcr, gcz, gcl
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: n0 = 0, m0 = 0, m1 = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsfix, iflag, js
      REAL(rprec) :: r1, r2, tnorm, fac, bz1
C-----------------------------------------------

!
!     IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!     INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!     ZCS = RSS, ZCC = RSC ARE THE CORRECT POLAR RELATIONS; HERE WE USE
!     THE SIMPLER RELATION ZCC = 0
!

      IF (lthreed .and. (fsqz.lt.1.e-8 .or. lprec2d))
     1   gcz(:,1:ntor,m1,zcs) = 0                                   !ZCS = 0
      IF (lasym) THEN
         gcr(:,:,m1,rsc) = p5*(gcr(:,:,m1,rsc) + gcz(:,:,m1,zcc))
         gcz(:,:,m1,zcc) = gcr(:,:,m1,rsc)                          !ZCC = RSC
         IF (lprec2d) gcz(:,:,m1,zcc) = 0                           !Avoid singular Hessian
         IF (lthreed) THEN
            gcr(:,:,m1,rcs) = p5*(gcr(:,:,m1,rcs) - gcz(:,:,m1,zss))
            gcz(:,:,m1,zss) = gcr(:,:,m1,rcs)
            !Not implemented yet in Hessian; set gcz(zss) = 0 and tie blocks as for zcc above, with - sign
         END IF
      END IF

      IF (iequi .ge. 2) RETURN

!FOR DEBUGGING:      IF (lprec2d) gcz(:,0,0,zcc) = 0
!
!     PUT FORCES INTO PHIFSAVE UNITS THAT PRECONDITIONERS, FNORM ARE IN
!
      IF (phifac .eq. zero) THEN
         STOP 'phifac = 0 in residue'
      ELSE
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      END IF


      IF (lrecon) THEN
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS

        r1 = SUM(gcr(:ns,n0,m0,1))
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        IF (iopt_raxis.gt.0 .and. iresidue.eq.2
     1     .and. fsq.lt.fopt_axis) iresidue = 3
        IF (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      ELSE
!
!     ADJUST PHIEDGE
!
         IF (imovephi .gt. 0) CALL movephi1 (gphifac)
      ENDIF
      gc(neqs1) = gphifac

!
!       CONSTRUCT INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, m0)
      bz1 = 2*hs*r1/(rbtor*tnorm*phip(2))**2
      fsql = bz1*SUM(gcl*gcl)

      r2 = SUM(gcr(ns,:,:,:)*gcr(ns,:,:,:) 
     1   +     gcz(ns,:,:,:)*gcz(ns,:,:,:))
      fedge = r1*fnorm*r2
!
!       PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      IF (.not.lprec2d) THEN
         iflag = 0
         CALL scalfor (gcr, arm, brm, ard, brd, crd, iflag)
         iflag = 1
         CALL scalfor (gcz, azm, bzm, azd, bzd, crd, iflag)

         gcl = faclam*gcl*r1

         IF (lasym) THEN
            DO js = 2, ns
               gcz(js,:,m1,zcc) = (ard(js,2)*gcr(js,:,m1,rsc)
     1                          +  azd(js,2)*gcz(js,:,m1,zcc))
               gcz(js,:,m1,zcc) = gcz(js,:,m1,zcc)/(ard(js,2)+azd(js,2))
               gcr(js,:,m1,rsc) = gcz(js,:,m1,zcc)
            END DO
         END IF

         CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1/4, m1)
 
         fac = .5_dp
         IF ((iter2-iter1).lt.25 .and. (fsqr+fsqz).gt.1.E-2_dp)
     1     fac = fac / sqrt(1.E2_dp*(fsqr+fsqz))
         gcr = fac*gcr
         gcz = fac*gcz

!
!     SOFT-STARTUP AVOIDS SOME "HANG-UPS"  (04/02)
!
         IF (iter2.lt.50 .and. (fsqr+fsqz).gt.1.E-3_dp) 
     1       gcl = gcl/3
         fsql1 = hs*SUM(gcl*gcl)/r1

      ELSE
         CALL block_precond(gc)
         IF (.not.lfreeb .and. ANY(gcr(ns,:,:,:) .ne. zero))
     1      STOP 'gcr(ns) != 0 for fixed boundary in residue'
         IF (.not.lfreeb .and. ANY(gcz(ns,:,:,:) .ne. zero))
     1      STOP 'gcz(ns) != 0 for fixed boundary in residue'
         IF (ANY(gcl(:,:,0,1) .ne. zero))
     1      STOP 'gcl(m=0) != 0'
         IF (lthreed) THEN
            IF (ANY(gcl(:,0,:,2) .ne. zero))
     1      STOP 'gcl(n=0) != 0'
         END IF
         IF (lasym) THEN
            IF (ANY(gcz(:,:,m1,zcc) .ne. zero))
     1      STOP 'gcz(m=1,cc) ! = 0 in residue!'
            gcz(:,:,m1,zcc) = gcr(:,:,m1,rsc)                        !ZCC(m=1) = RSC(m=1)
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ENDIF

      END SUBROUTINE residue
