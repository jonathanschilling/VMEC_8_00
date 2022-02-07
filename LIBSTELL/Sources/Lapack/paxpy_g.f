      SUBROUTINE PAXPY_G (N, SA, SX, SY, INDX)
      USE LAPREC
      IMPLICIT NONE
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION AX + Y(INDX(I) ROUTINES
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INDX(N)
      REAL(RPREC) :: SA, SX(N), SY(*)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: ZERO = 0
C-----------------------------------------------

      IF (n .le. 0) RETURN
      IF (sa .eq. zero) RETURN

      sy(indx(:n)) = sy(indx(:n)) + sa * sx(:n)

      END SUBROUTINE PAXPY_G
