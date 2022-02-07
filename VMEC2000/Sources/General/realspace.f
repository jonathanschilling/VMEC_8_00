      MODULE realspace
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   r1, ru, rv, zu, zv, rcon, zcon
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu, guv,
     1   gvv, ru0, zu0, rcon0, zcon0, phip, shalf, sqrts, wint
C-----------------------------------------------
      END MODULE realspace
