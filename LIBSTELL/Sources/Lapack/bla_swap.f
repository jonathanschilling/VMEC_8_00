c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE sswap (N, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), INTENT(INOUT) :: SX(*), SY(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
      REAL(WP), DIMENSION(:), ALLOCATABLE :: STEMP
C-----------------------------------------------
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      IF (n .le. 0) RETURN

      ALLOCATE (stemp(n), stat=ix)
      IF (ix .ne. 0) STOP 'Allocation error in swap'

      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c       code for unequal increments or equal increments not equal
c         to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         IF (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         stemp(:n) = sx(ix:(n-1)*incx+ix:incx)
         sx(ix:(n-1)*incx+ix:incx) = sy(iy:(n-1)*incy+iy:incy)
         sy(iy:(n-1)*incy+iy:incy) = stemp(:n)

      ELSE
c
c     code for both increments equal to 1
c
         stemp(:n) = sx(:n)
         sx(:n) = sy(:n)
         sy(:n) = stemp(:n)
      END IF

      DEALLOCATE (stemp)

      END SUBROUTINE sswap

      SUBROUTINE dswap (N, SX, INCX, SY, INCY)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), INTENT(INOUT) :: SX(*), SY(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
      REAL(WP), DIMENSION(:), ALLOCATABLE :: STEMP
C-----------------------------------------------
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      IF (n .le. 0) RETURN

      ALLOCATE (stemp(n), stat=ix)
      IF (ix .ne. 0) STOP 'Allocation error in swap'

      IF (incx.ne.1 .or. incy.ne.1) THEN
c
c       code for unequal increments or equal increments not equal
c         to 1
c
         ix = 1
         iy = 1
         IF (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         IF (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         stemp(:n) = sx(ix:(n-1)*incx+ix:incx)
         sx(ix:(n-1)*incx+ix:incx) = sy(iy:(n-1)*incy+iy:incy)
         sy(iy:(n-1)*incy+iy:incy) = stemp(:n)

      ELSE
c
c     code for both increments equal to 1
c
         stemp(:n) = sx(:n)
         sx(:n) = sy(:n)
         sy(:n) = stemp(:n)
      END IF

      DEALLOCATE (stemp)

      END SUBROUTINE dswap
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
