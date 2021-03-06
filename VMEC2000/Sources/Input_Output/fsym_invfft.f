      SUBROUTINE fsym_invfft (bsubsu, bsubsv, ns)
      USE vmec_main, ONLY: rprec, nzeta, ntheta1, ntheta2, 
     1    ntheta3, ireflect
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, jkz, jkr
C-----------------------------------------------

      DO i = 1 + ntheta2, ntheta1
         ir = ntheta1 + 2 - i                 !-theta
         DO jkz= 1, ns*nzeta
            jkr = ireflect(jkz)               !-zeta
            bsubsu(jkz,i,0) = bsubsu(jkr,ir,0) - bsubsu(jkr,ir,1)
            bsubsv(jkz,i,0) = bsubsv(jkr,ir,0) - bsubsv(jkr,ir,1)
         END DO
      END DO

      bsubsu(:,:ntheta2,0)=bsubsu(:,:ntheta2,0) + bsubsu(:,:ntheta2,1)
      bsubsv(:,:ntheta2,0)=bsubsv(:,:ntheta2,0) + bsubsv(:,:ntheta2,1)

      END SUBROUTINE fsym_invfft
