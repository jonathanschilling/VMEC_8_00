c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE sscal (N, SA, SX, INCX)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX
      REAL(WP), INTENT(IN) :: SA
      REAL(WP), INTENT(INOUT) :: SX(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX
C-----------------------------------------------
!
!     scales a vector by a constant.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increments, 9/29/88.
!
      IF (n .le. 0) RETURN
      IF (incx .ne. 1) THEN
!
!     code for unequal increments or equal increments
!     not equal to 1
!
         ix = 1
         IF(incx.lt.0)ix = (-n+1)*incx + 1
         sx(ix:(n-1)*incx+ix:incx) = sa*sx(ix:(n-1)*incx+ix:incx)
      ELSE
!
!     code for both increments equal to 1
!
         sx(:n) = sa * sx(:n)
      END IF

      END SUBROUTINE sscal

      SUBROUTINE dscal (N, SA, SX, INCX)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX
      REAL(WP), INTENT(IN) :: SA
      REAL(WP), INTENT(INOUT) :: SX(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX
C-----------------------------------------------
!
!     scales a vector by a constant.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increments, 9/29/88.
!
      IF (n .le. 0) RETURN
      IF (incx .ne. 1) THEN
!
!     code for unequal increments or equal increments
!     not equal to 1
!
         ix = 1
         IF(incx.lt.0)ix = (-n+1)*incx + 1
         sx(ix:(n-1)*incx+ix:incx) = sa*sx(ix:(n-1)*incx+ix:incx)
      ELSE
!
!     code for both increments equal to 1
!
         sx(:n) = sa * sx(:n)
      END IF

      END SUBROUTINE dscal
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
