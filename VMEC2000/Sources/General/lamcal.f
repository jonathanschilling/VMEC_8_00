      SUBROUTINE lamcal(phipog, guu, guv, gvv)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, jlam
      USE realspace, ONLY: sqrts, phip
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1   phipog, guu, guv, gvv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: damping_fac = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m,n,js
      REAL(rprec) :: tnn, tnm, tmm, power
C-----------------------------------------------

      blam(:ns) = SUM(guu*phipog, dim=2)*phip(2)
      clam(:ns) = SUM(gvv*phipog, dim=2)*phip(2)
      dlam(:ns) = SUM(guv*phipog, dim=2)*phip(2)
      blam(1) =  blam(2)
      clam(1) =  clam(2)
      dlam(1) =  dlam(2)
      blam(ns+1) =  0
      clam(ns+1) =  0
      dlam(ns+1) =  0
      DO js = 2, ns
        blam(js) = cp5*(blam(js) + blam(js+1))
        clam(js) = cp5*(clam(js) + clam(js+1))
        dlam(js) = cp5*(dlam(js) + dlam(js+1))
      END DO

      faclam = 0
      DO m = 0, mpol1
         tmm = m*m
         power = MIN(tmm/256, 8._dp)
         DO n = 0, ntor
            IF (m.eq.0 .and. n.eq.0) CYCLE
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            DO js = jlam(m), ns
               faclam(js,n,m,1) = 2*damping_fac/
     1               ((blam(js) + blam(js+1))*tnn
     2         + SIGN((dlam(js) + dlam(js+1)),blam(js))*tnm
     3         +      (clam(js) + clam(js+1))*tmm)
     4         * sqrts(js)**power                                   !Damps m > 16 modes
            END DO
         END DO
      END DO

      DO n = 2, ntmax
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
      END DO

!
!     ADD NORM FOR IOTA FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
!
      IF (liota) THEN
         DO js = 1, ns
            faclam(js,0,0,1) = 2*damping_fac/(blam(js) + blam(js+1))
         END DO
      END IF
      END SUBROUTINE lamcal
