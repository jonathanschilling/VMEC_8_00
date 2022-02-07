      SUBROUTINE splINT(xa, ya, y2a, n, x, y, yp, ndim)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, ndim
      REAL(rprec), DIMENSION(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: klo, khi, i, k
      REAL(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, deriv
C-----------------------------------------------
!
!       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
!       XA: ordered array of length N of ordinates at which FUNCTION YA=F(XA)
!           is tabulated
!       YA: array of length N , = F(XA)
!       Y2A: array of second derivatives at XA points
!       computed from CALL to SPLINE
!       X : value at which Y = F(X) is to be computed from splines
!       YP = dY/dX at X
!       NDIM: DIMENSION of X, Y, YP arrays


      deriv = yp(1)
      klo = 1
      khi = n
      DO i = 1, ndim
         DO WHILE(khi - klo .gt. 1)
            k = (khi + klo)/2
            IF (xa(k) .gt. x(i)) THEN
               khi = k
            ELSE
               klo = k
            ENDIF
         END DO

         h = xa(khi) - xa(klo)
         a = xa(khi) - x(i)
         b = x(i) - xa(klo)
         h2 = h*h
         a2 = a*a
         b2 = b*b
         y26lo = c1o6*y2a(klo)
         y26hi = c1o6*y2a(khi)
         y(i) = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h
         IF (deriv .ne. zero) yp(i) = (ya(khi)-ya(klo)+y26hi*(3.0*b2-h2)
     1      -y26lo*(3.0*a2-h2))/h
         IF (i.lt.ndim .and. x(i+1).gt.x(i)) THEN
            khi = n
         ELSE
            klo = 1
         ENDIF
      END DO

      END SUBROUTINE splint
