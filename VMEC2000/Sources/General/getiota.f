      SUBROUTINE getiota(phipog, bsupu, bsupv)
      USE vmec_main
      USE realspace
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in) :: phipog, bsupv
      REAL(rprec), DIMENSION(nrzt), INTENT(inout) :: bsupu
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l
      REAL(rprec) :: top, bot
C-----------------------------------------------

      IF (ncurr .eq. 0) RETURN
      DO js = 2, ns
         top = icurv(js)                                !!Integrated toroidal current
         bot = 0
         DO l = js, nrzt, ns
            top = top - wint(l)*(guu(l)*bsupu(l) + guv(l)*bsupv(l))
            bot = bot + wint(l)*phipog(l)*guu(l)
         END DO
         iotas(js) = top/bot
      END DO

!     Do not compute iota too near origin
      iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)           !!zero gradient near axis
!      iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
!      DO js = 2, ns-1
!         iotaf(js) = p5*(iotas(js) + iotas(js+1))
!      END DO
      DO js = 2, ns
         iotaf(js) = 2*iotas(js) - iotaf(js-1)
      END DO

      END SUBROUTINE getiota
