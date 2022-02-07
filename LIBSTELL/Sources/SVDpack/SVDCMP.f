      SUBROUTINE SVDCMP(A, M, N, MP, NP, W, V)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP) :: A
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(NP,NP) :: V
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0, one = 1, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, L, K, J, ITS, NM
      REAL(rprec), DIMENSION(MP*NP) :: RV1
      REAL(rprec) :: G, SCALE, ANORM, S, F, H, C, Y, Z, X
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PYTHAG
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SVDCMP'
      IF (N .GT. NP) STOP 'N > NP IN SVDCMP'

      G = zero
      SCALE = zero
      ANORM = zero
      DO I = 1, N
         L = I + 1
         RV1(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I .le. M) THEN
            SCALE = SUM(ABS(A(I:M,I)))
            IF (SCALE .ne. zero) THEN
               A(I:M,I) = A(I:M,I)/SCALE
               S = SUM(A(I:M,I)*A(I:M,I))
               F = A(I,I)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A(I,I) = F - G
               IF (I .ne. N) THEN
                  DO J = L, N
                     S = SUM(A(I:M,I)*A(I:M,J))
                     F = S/H
                     A(I:M,J) = A(I:M,J) + F*A(I:M,I)
                  END DO
               ENDIF
               A(I:M,I) = SCALE*A(I:M,I)
            ENDIF
         ENDIF
         W(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I.le.M .AND. I.ne.N) THEN
            SCALE = SUM(ABS(A(I,L:N)))
            IF (SCALE .ne. zero) THEN
               A(I,L:N) = A(I,L:N)/SCALE
               S = SUM(A(I,L:N)*A(I,L:N))
               F = A(I,L)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A(I,L) = F - G
               RV1(L:N) = A(I,L:N)/H
               IF (I .ne. M) THEN
                  DO J = L, M
                     S = SUM(A(J,L:N)*A(I,L:N))
                     A(J,L:N) = A(J,L:N) + S*RV1(L:N)
                  END DO
               ENDIF
               A(I,L:N) = SCALE*A(I,L:N)
            ENDIF
         ENDIF
         ANORM = MAX(ANORM,ABS(W(I))+ABS(RV1(I)))
      END DO
      DO I = N, 1, -1
         IF (I .lt. N) THEN
            IF (G .ne. zero) THEN
               V(L:N,I) = (A(I,L:N)/A(I,L))/G
               DO J = L, N
                  S = SUM(A(I,L:N)*V(L:N,J))
                  V(L:N,J) = V(L:N,J) + S*V(L:N,I)
               END DO
            ENDIF
            V(I,L:N) = zero
            V(L:N,I) = zero
         ENDIF
         V(I,I) = one
         G = RV1(I)
         L = I
      END DO
      DO I = N, 1, -1
         L = I + 1
         G = W(I)
         IF (I .lt. N) THEN
            A(I,L:N) = zero
         ENDIF
         IF (G .ne. zero) THEN
            G = one/G
            IF (I .ne. N) THEN
               DO J = L, N
                  S = SUM(A(L:M,I)*A(L:M,J))
                  F = (S/A(I,I))*G
                  A(I:M,J) = A(I:M,J) + F*A(I:M,I)
               END DO
            ENDIF
            A(I:M,I) = A(I:M,I)*G
         ELSE
            A(I:M,I) = zero
         ENDIF
         A(I,I) = A(I,I) + one
      END DO
      L49: DO K = N, 1, -1
         DO ITS = 1, 50
            DO L = K, 1, -1
               NM = L - 1
               IF (ABS(RV1(L)) + ANORM .eq. ANORM) GO TO 2
               IF (ABS(W(NM)) + ANORM .eq. ANORM) EXIT
            END DO
            C = zero
            S = one
            DO I = L, K
               F = S*RV1(I)
               IF (ABS(F) + ANORM .ne. ANORM) THEN
                  G = W(I)
                  H = PYTHAG(F,G)
c              H=SQRT(F*F+G*G)
                  W(I) = H
                  H = one/H
                  C = G*H
                  S = -F*H
                  DO J = 1, M
                     Y = A(J,NM)
                     Z = A(J,I)
                     A(J,NM) = Y*C + Z*S
                     A(J,I) = (-Y*S) + Z*C
                  END DO
               ENDIF
            END DO
    2       CONTINUE
            Z = W(K)
            IF (L .eq. K) THEN
               IF (Z .lt. zero) THEN
                  W(K) = -Z
                  V(:N,K) = -V(:N,K)
               ENDIF
               CYCLE  L49
            ENDIF
            IF (ITS .eq. 50) THEN
               WRITE (*, '(2A)') 'PAUSE ',
     1            'No convergence in 150 iterations'
               READ *
            ENDIF
            X = W(L)
            NM = K - 1
            Y = W(NM)
            G = RV1(NM)
            H = RV1(K)
            F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(TWO*H*Y)
            G = PYTHAG(F,one)
c          G=SQRT(F*F + one)
            F = ((X - Z)*(X + Z) + H*(Y/(F + SIGN(G,F)) - H))/X
            C = one
            S = one
            DO J = L, NM
               I = J + 1
               G = RV1(I)
               Y = W(I)
               H = S*G
               G = C*G
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               RV1(J) = Z
               C = F/Z
               S = H/Z
               F = X*C + G*S
               G = (-X*S) + G*C
               H = Y*S
               Y = Y*C
               DO NM = 1, N
                  X = V(NM,J)
                  Z = V(NM,I)
                  V(NM,J) = X*C + Z*S
                  V(NM,I) = (-X*S) + Z*C
               END DO
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               W(J) = Z
               IF (Z .ne. zero) THEN
                  Z = one/Z
                  C = F*Z
                  S = H*Z
               ENDIF
               F = C*G + S*Y
               X = (-S*G) + C*Y
               DO NM = 1, M
                  Y = A(NM,J)
                  Z = A(NM,I)
                  A(NM,J) = Y*C + Z*S
                  A(NM,I) = (-Y*S) + Z*C
               END DO
            END DO
            RV1(L) = zero
            RV1(K) = F
            W(K) = X
         END DO
      END DO L49

      END SUBROUTINE SVDCMP
