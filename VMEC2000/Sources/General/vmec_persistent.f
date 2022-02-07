      MODULE vmec_persistent
      USE stel_kinds, ONLY: rprec1 => rprec
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixm, jmin3
      REAL(rprec1), DIMENSION(:,:), ALLOCATABLE :: cosmu, sinmu,
     1   cosmum, sinmum, cosmui, cosmumi, sinmui, sinmumi,
     2   cosnv, sinnv, cosnvn, sinnvn
      REAL(rprec1), DIMENSION(:), ALLOCATABLE ::
     1      xm, xn, xm_nyq, xn_nyq, cos01, sin01
c-----------------------------------------------
      END MODULE vmec_persistent
