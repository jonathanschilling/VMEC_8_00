      SUBROUTINE SVBKSB(U, W, V, M, N, MP, NP, B, X)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP) :: U
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(NP,NP) :: V
      REAL(rprec), DIMENSION(MP) :: B
      REAL(rprec), DIMENSION(NP) :: X
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J
      REAL(rprec), DIMENSION(MP*NP) :: TMP
      REAL(rprec) :: S
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SVBKSB'
      IF (N .GT. NP) STOP 'N > NP IN SVBKSB'

      DO J = 1, N
         S = zero
         IF (W(J) .ne. zero) THEN
            S = SUM(U(:M,J)*B(:M))
            S = S/W(J)
         ENDIF
         TMP(J) = S
      END DO
      DO J = 1, N
         S = SUM(V(J,:N)*TMP(:N))
         X(J) = S
      END DO

      END SUBROUTINE SVBKSB
