      SUBROUTINE solver(amat, b, m)
      USE stel_kinds
      USE laprec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m
      REAL(rprec), INTENT(inout) :: amat(m,*), b(m,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: nrhs = 1
      INTEGER :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
C-----------------------------------------------
      info = 0
      ALLOCATE (ipiv(m))

!     Compute the solution to a REAL system of linear equations
!       A * X = B,
!     WHERE A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!
!     FACTOR AMATRIX INTO LU FORM
!     AND SOLVE BY GAUSSIAN ELIMINATION
!
      CALL la_gesv (m, nrhs, amat, m, ipiv, b, m, info)
      IF (info .ne. 0) PRINT *, ' Condition No. = 0 in Solver'

      DEALLOCATE (ipiv)

      END SUBROUTINE solver
