c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE strsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1    A, LDA, B, LDB)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, LDA, LDB
      REAL(WP), INTENT(IN) :: ALPHA
      CHARACTER, INTENT(IN) :: SIDE, UPLO, TRANSA, DIAG
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB,*), INTENT(INOUT) :: B
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: I, INFO, J, K, NROWA
      REAL(WP) :: TEMP
      LOGICAL :: LSIDE, NOUNIT, UPPER
C-----------------------------------------------
C   I N T R I N S I C  F U N C T I O N S
C-----------------------------------------------
      INTRINSIC MAX
      LOGICAL, EXTERNAL :: lsame
C-----------------------------------------------

*     .. ARRAY ARGUMENTS ..
*     ..
*
*  Purpose
*  =======
*
*  STRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  WHERE alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A''.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On ENTRY, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on EXIT.
*
*  UPLO   - CHARACTER*1.
*           On ENTRY, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on EXIT.
*
*  TRANSA - CHARACTER*1.
*           On ENTRY, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A''.
*
*              TRANSA = 'C' or 'c'   op( A ) = A''.
*
*           Unchanged on EXIT.
*
*  DIAG   - CHARACTER*1.
*           On ENTRY, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on EXIT.
*
*  M      - INTEGER.
*           On ENTRY, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero THEN  A is not referenced and  B need not be set before
*           ENTRY.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, k ), WHERE k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before ENTRY  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before ENTRY  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on EXIT.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  THEN
*           LDA  must be at least  MAX( 1, m ),  when  SIDE = 'R' or 'r'
*           THEN LDA must be at least MAX( 1, n ).
*           Unchanged on EXIT.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before ENTRY,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on EXIT  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On ENTRY, LDB specifies the first DIMENSION of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           MAX( 1, m ).
*           Unchanged on EXIT.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
*     .. Local Scalars ..
*     .. Parameters ..
      LSIDE = lsame(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      ENDIF
      NOUNIT = lsame(DIAG,'N')
      UPPER = lsame(UPLO,'U')
*
      INFO = 0
      IF ( .NOT.LSIDE .AND.  .NOT.lsame(SIDE,'R')) THEN
         INFO = 1
      ELSE IF ( .NOT.UPPER .AND.  .NOT.lsame(UPLO,'L')) THEN
         INFO = 2
      ELSE IF (.NOT.lsame(TRANSA,'N') .AND. .NOT.lsame(TRANSA,'T')
     1       .AND.  .NOT.lsame(TRANSA,'C')) THEN
         INFO = 3
      ELSE IF ( .NOT.lsame(DIAG,'U') .AND.
     1          .NOT.lsame(DIAG,'N')) THEN
         INFO = 4
      ELSE IF (M < 0) THEN
         INFO = 5
      ELSE IF (N < 0) THEN
         INFO = 6
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB < MAX(1,M)) THEN
         INFO = 11
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN STRSM '
         RETURN
      ENDIF
*
*     Quick RETURN IF possible.
*
      IF (N .EQ. 0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         B(:M,:N) = ZERO
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (LSIDE) THEN
         IF (lsame(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = M, 1, -1
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(:K-1,J) = B(:K-1,J) - B(K,J)*A(:K-1,K)
                     ENDIF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, M
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(K+1:M,J) = B(K+1:M,J) - B(K,J)*A(K+1:M,K)
                     ENDIF
                  END DO
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*inv( A'' )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  DO I = 1, M
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(:I-1,I)*B(:I-1,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = M, 1, -1
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(I+1:M,I)*B(I+1:M,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ENDIF
         ENDIF
      ELSE
         IF (lsame(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, J - 1
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ELSE
               DO J = N, 1, -1
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = J + 1, N
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*B*inv( A'' ).
*
            IF (UPPER) THEN
               DO K = N, 1, -1
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = 1, K - 1
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ELSE
               DO K = 1, N
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = K + 1, N
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE strsm

      SUBROUTINE dtrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1    A, LDA, B, LDB)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, LDA, LDB
      REAL(WP), INTENT(IN) :: ALPHA
      CHARACTER, INTENT(IN) :: SIDE, UPLO, TRANSA, DIAG
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB,*), INTENT(INOUT) :: B
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: I, INFO, J, K, NROWA
      REAL(WP) :: TEMP
      LOGICAL :: LSIDE, NOUNIT, UPPER
C-----------------------------------------------
C   I N T R I N S I C  F U N C T I O N S
C-----------------------------------------------
      INTRINSIC MAX
      LOGICAL, EXTERNAL :: lsame
C-----------------------------------------------

*     .. ARRAY ARGUMENTS ..
*     ..
*
*  Purpose
*  =======
*
*  STRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  WHERE alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A''.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On ENTRY, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on EXIT.
*
*  UPLO   - CHARACTER*1.
*           On ENTRY, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on EXIT.
*
*  TRANSA - CHARACTER*1.
*           On ENTRY, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A''.
*
*              TRANSA = 'C' or 'c'   op( A ) = A''.
*
*           Unchanged on EXIT.
*
*  DIAG   - CHARACTER*1.
*           On ENTRY, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on EXIT.
*
*  M      - INTEGER.
*           On ENTRY, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero THEN  A is not referenced and  B need not be set before
*           ENTRY.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, k ), WHERE k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before ENTRY  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before ENTRY  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on EXIT.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  THEN
*           LDA  must be at least  MAX( 1, m ),  when  SIDE = 'R' or 'r'
*           THEN LDA must be at least MAX( 1, n ).
*           Unchanged on EXIT.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before ENTRY,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on EXIT  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On ENTRY, LDB specifies the first DIMENSION of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           MAX( 1, m ).
*           Unchanged on EXIT.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
*     .. Local Scalars ..
*     .. Parameters ..
      LSIDE = lsame(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      ENDIF
      NOUNIT = lsame(DIAG,'N')
      UPPER = lsame(UPLO,'U')
*
      INFO = 0
      IF ( .NOT.LSIDE .AND.  .NOT.lsame(SIDE,'R')) THEN
         INFO = 1
      ELSE IF ( .NOT.UPPER .AND.  .NOT.lsame(UPLO,'L')) THEN
         INFO = 2
      ELSE IF (.NOT.lsame(TRANSA,'N') .AND. .NOT.lsame(TRANSA,'T')
     1       .AND.  .NOT.lsame(TRANSA,'C')) THEN
         INFO = 3
      ELSE IF ( .NOT.lsame(DIAG,'U') .AND.
     1          .NOT.lsame(DIAG,'N')) THEN
         INFO = 4
      ELSE IF (M < 0) THEN
         INFO = 5
      ELSE IF (N < 0) THEN
         INFO = 6
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB < MAX(1,M)) THEN
         INFO = 11
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN STRSM '
         RETURN
      ENDIF
*
*     Quick return if possible.
*
      IF (N .EQ. 0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         B(:M,:N) = ZERO
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (LSIDE) THEN
         IF (lsame(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = M, 1, -1
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(:K-1,J) = B(:K-1,J) - B(K,J)*A(:K-1,K)
                     ENDIF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, M
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(K+1:M,J) = B(K+1:M,J) - B(K,J)*A(K+1:M,K)
                     ENDIF
                  END DO
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*inv( A'' )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  DO I = 1, M
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(:I-1,I)*B(:I-1,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = M, 1, -1
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(I+1:M,I)*B(I+1:M,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ENDIF
         ENDIF
      ELSE
         IF (lsame(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, J - 1
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ELSE
               DO J = N, 1, -1
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = J + 1, N
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*B*inv( A'' ).
*
            IF (UPPER) THEN
               DO K = N, 1, -1
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = 1, K - 1
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ELSE
               DO K = 1, N
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = K + 1, N
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE dtrsm
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
