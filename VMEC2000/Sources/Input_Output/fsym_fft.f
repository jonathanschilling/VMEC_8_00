      SUBROUTINE fsym_fft (bs, bu, bv, bs_s, bu_s, bv_s,
     1    bs_a, bu_a, bv_a)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(in) :: bs
      REAL(rprec), DIMENSION(nzeta,ntheta3,0:1), INTENT(in) :: bu, bv
      REAL(rprec), DIMENSION(nzeta,ntheta2,0:1), INTENT(out) ::
     1   bs_s, bu_s, bv_s, bs_a, bu_a, bv_a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, mpar, kz, kr
C-----------------------------------------------

!
!       SYMMETRIZE CONTRAVARIANT COMPONENTS OF B
!       SO COS,SIN INTEGRALS CAN BE PERFORMED ON HALF-THETA INTERVAL
!
!       bs_s(v,u) = .5*( bs_s(v,u) - bs_s(-v,-u) )     ! * SIN(mu - nv)
!       bs_a(v,u) = .5*( bs_s(v,u) + bs_s(-v,-u) )     ! * COS(mu - nv)
!
!
      DO mpar = 0, 1
         DO i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            IF (i == 1) ir = 1
            DO kz = 1, nzeta
               kr = ireflect(ns*kz)/ns           !-zeta
               bs_a(kz,i,mpar) = cp5*(bs(kz,i)+bs(kr,ir))
               bs_s(kz,i,mpar) = cp5*(bs(kz,i)-bs(kr,ir))
               bu_a(kz,i,mpar) = cp5*(bu(kz,i,mpar)-bu(kr,ir,mpar))
               bu_s(kz,i,mpar) = cp5*(bu(kz,i,mpar)+bu(kr,ir,mpar))
               bv_a(kz,i,mpar) = cp5*(bv(kz,i,mpar)-bv(kr,ir,mpar))
               bv_s(kz,i,mpar) = cp5*(bv(kz,i,mpar)+bv(kr,ir,mpar))
            END DO
         END DO
      END DO

      END SUBROUTINE fsym_fft
