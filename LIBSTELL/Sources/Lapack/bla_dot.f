c !DEC$ IF DEFINED (NEED_BLAS)
      FUNCTION sdot (n, sx, incx, sy, incy)
      USE laprec, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: n, incx, incy
      REAL(WP), INTENT(IN) :: sx(*), sy(*)
      REAL(WP) :: sdot
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
c
      IF (n .le. 0) RETURN
      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = (-n+1)*incx + 1
         IF (incy .lt. 0) iy = (-n+1)*incy + 1
         sdot = SUM (sx(ix:(n-1)*incx+ix:incx)*
     1                   sy(iy:(n-1)*incy+iy:incy))
      ELSE
c
c     code for both increments equal to 1
c
         sdot = SUM(sx(:n)*sy(:n))
      END IF

      END FUNCTION sdot

      FUNCTION ddot (n, sx, incx, sy, incy)
      USE laprec, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: n, incx, incy
      REAL(WP), INTENT(IN) :: sx(*), sy(*)
      REAL(WP) :: ddot
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
c
      IF (n .le. 0) RETURN
      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = (-n+1)*incx + 1
         IF (incy .lt. 0) iy = (-n+1)*incy + 1
         ddot = SUM (sx(ix:(n-1)*incx+ix:incx)*
     1                   sy(iy:(n-1)*incy+iy:incy))
      ELSE
c
c     code for both increments equal to 1
c
         ddot = SUM(sx(:n)*sy(:n))
      END IF

      END FUNCTION ddot
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
