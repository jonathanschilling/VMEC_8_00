c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE saxpy(N, SA, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), INTENT(IN) :: SA, SX(*)
      REAL(WP), INTENT(INOUT) :: SY(*)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
C
      IF (n .le. 0) RETURN
      IF (sa .eq. zero) RETURN

      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0)ix = (-n+1)*incx + 1
         IF (incy .lt. 0)iy = (-n+1)*incy + 1

         sy(iy:(n-1)*incy+iy:incy) = sy(iy:(n-1)*incy+iy:incy)
     1            + sa * sx(ix:(n-1)*incx+ix:incx)
      ELSE
c
c     code for both increments equal to 1
c
         sy(:n) = sy(:n) + sa*sx(:n)
      END IF

      END SUBROUTINE saxpy

      SUBROUTINE daxpy(N, SA, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), INTENT(IN) :: SA, SX(*)
      REAL(WP), INTENT(INOUT) :: SY(*)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
C
      IF (n .le. 0) RETURN
      IF (sa .eq. zero) RETURN

      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0)ix = (-n+1)*incx + 1
         IF (incy .lt. 0)iy = (-n+1)*incy + 1

         sy(iy:(n-1)*incy+iy:incy) = sy(iy:(n-1)*incy+iy:incy)
     1            + sa * sx(ix:(n-1)*incx+ix:incx)
      ELSE
c
c     code for both increments equal to 1
c
         sy(:n) = sy(:n) + sa*sx(:n)
      END IF

      END SUBROUTINE daxpy
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
