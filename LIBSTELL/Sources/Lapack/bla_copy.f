c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE scopy (N, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), DIMENSION(*), INTENT(IN) :: SX
      REAL(WP), DIMENSION(*), INTENT(OUT) :: SY
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      IF (n .le. 0) RETURN
      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c        code for unequal increments or equal increments
c          not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         IF (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         sy(iy:(n-1)*incy+iy:incy) = sx(ix:(n-1)*incx+ix:incx)
      ELSE
c
c     code for both increments equal to 1
c
         sy(:n) = sx(:n)
      END IF

      END SUBROUTINE scopy

      SUBROUTINE dcopy (N, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), DIMENSION(*), INTENT(IN) :: SX
      REAL(WP), DIMENSION(*), INTENT(OUT) :: SY
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      IF (n .le. 0) RETURN
      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c        code for unequal increments or equal increments
c          not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         IF (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         sy(iy:(n-1)*incy+iy:incy) = sx(ix:(n-1)*incx+ix:incx)
      ELSE
c
c     code for both increments equal to 1
c
         sy(:n) = sx(:n)
      END IF

      END SUBROUTINE dcopy
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
