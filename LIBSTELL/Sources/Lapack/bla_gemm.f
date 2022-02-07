c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE sgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     1                  B, LDB, BETA, C, LDC )
      USE LAPREC, ONLY: WP => SP
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      CHARACTER*1         TRANSA, TRANSB
      INTEGER             M, N, K, LDA, LDB, LDC
      REAL(WP) :: ALPHA, BETA
      REAL(WP) :: A( LDA, * ), B( LDB, * ), C( LDC, * )
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(WP) :: TEMP
      LOGICAL :: NOTA, NOTB
      LOGICAL, EXTERNAL :: lsame
C-----------------------------------------------
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  WHERE  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On ENTRY, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A''.
*
*              TRANSA = 'C' or 'c',  op( A ) = A''.
*
*           Unchanged on EXIT.
*
*  TRANSB - CHARACTER*1.
*           On ENTRY, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B''.
*
*              TRANSB = 'C' or 'c',  op( B ) = B''.
*
*           Unchanged on EXIT.
*
*  M      - INTEGER.
*           On ENTRY,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on EXIT.
*
*  K      - INTEGER.
*           On ENTRY,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY, ALPHA specifies the scalar alpha.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), WHERE ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before ENTRY with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on EXIT.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) PROGRAM. When  TRANSA = 'N' or 'n' THEN
*           LDA must be at least  MAX( 1, m ), otherwise  LDA must be at
*           least  MAX( 1, k ).
*           Unchanged on EXIT.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), WHERE kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before ENTRY with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on EXIT.
*
*  LDB    - INTEGER.
*           On ENTRY, LDB specifies the first DIMENSION of B as declared
*           in the calling (sub) PROGRAM. When  TRANSB = 'N' or 'n' THEN
*           LDB must be at least  MAX( 1, k ), otherwise  LDB must be at
*           least  MAX( 1, n ).
*           Unchanged on EXIT.
*
*  BETA   - REAL            .
*           On ENTRY,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero THEN C need not be set on input.
*           Unchanged on EXIT.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before ENTRY, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           CASE C need not be set on ENTRY.
*           On EXIT, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On ENTRY, LDC specifies the first DIMENSION of C as declared
*           in  the  calling  (sub)  PROGRAM.   LDC  must  be  at  least
*           MAX( 1, m ).
*           Unchanged on EXIT.
*
*
*  Level 3 Blas routine.
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
      NOTA = lsame(TRANSA,'N')
      NOTB = lsame(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      ENDIF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      ENDIF
*
*     Test the input parameters.
*
      INFO = 0
      IF ( .NOT.NOTA .AND.  .NOT.lsame(TRANSA,'C') .AND. .
     1      NOT.lsame(TRANSA,'T')) THEN
         INFO = 1
      ELSE IF ( .NOT.NOTB .AND.  .NOT.lsame(TRANSB,'C') .AND.  .NOT.
     1      lsame(TRANSB,'T')) THEN
         INFO = 2
      ELSE IF (M < 0) THEN
         INFO = 3
      ELSE IF (N < 0) THEN
         INFO = 4
      ELSE IF (K < 0) THEN
         INFO = 5
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB < MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC < MAX(1,M)) THEN
         INFO = 13
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO, ' IN GEMM '
         RETURN
      ENDIF
*
*     Quick RETURN IF possible.
*
      IF(M.EQ.0.OR.N.EQ.0.OR.(ALPHA.EQ.ZERO.OR.K.EQ.0).AND.BETA.EQ.ONE)
     1   RETURN
*
*     And IF  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         IF (BETA .EQ. ZERO) THEN
            C(:M,:N) = ZERO
         ELSE
            C(:M,:N) = BETA*C(:M,:N)
         ENDIF
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (NOTB) THEN
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(L,J) .NE. ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(:K,J))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ELSE
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B'' + beta*C
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(J,L) .NE. ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B'' + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(J,:K))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF

      END SUBROUTINE sgemm

      SUBROUTINE dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     1                  B, LDB, BETA, C, LDC )
      USE LAPREC, ONLY: WP => DP
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      CHARACTER*1         TRANSA, TRANSB
      INTEGER             M, N, K, LDA, LDB, LDC
      REAL(WP) :: ALPHA, BETA
      REAL(WP) :: A( LDA, * ), B( LDB, * ), C( LDC, * )
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(WP) :: TEMP
      LOGICAL :: NOTA, NOTB
      LOGICAL, EXTERNAL :: lsame
C-----------------------------------------------
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  WHERE  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On ENTRY, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A''.
*
*              TRANSA = 'C' or 'c',  op( A ) = A''.
*
*           Unchanged on EXIT.
*
*  TRANSB - CHARACTER*1.
*           On ENTRY, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B''.
*
*              TRANSB = 'C' or 'c',  op( B ) = B''.
*
*           Unchanged on EXIT.
*
*  M      - INTEGER.
*           On ENTRY,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on EXIT.
*
*  K      - INTEGER.
*           On ENTRY,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY, ALPHA specifies the scalar alpha.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), WHERE ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before ENTRY with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on EXIT.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) PROGRAM. When  TRANSA = 'N' or 'n' THEN
*           LDA must be at least  MAX( 1, m ), otherwise  LDA must be at
*           least  MAX( 1, k ).
*           Unchanged on EXIT.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), WHERE kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before ENTRY with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on EXIT.
*
*  LDB    - INTEGER.
*           On ENTRY, LDB specifies the first DIMENSION of B as declared
*           in the calling (sub) PROGRAM. When  TRANSB = 'N' or 'n' THEN
*           LDB must be at least  MAX( 1, k ), otherwise  LDB must be at
*           least  MAX( 1, n ).
*           Unchanged on EXIT.
*
*  BETA   - REAL            .
*           On ENTRY,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero THEN C need not be set on input.
*           Unchanged on EXIT.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before ENTRY, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           CASE C need not be set on ENTRY.
*           On EXIT, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On ENTRY, LDC specifies the first DIMENSION of C as declared
*           in  the  calling  (sub)  PROGRAM.   LDC  must  be  at  least
*           MAX( 1, m ).
*           Unchanged on EXIT.
*
*
*  Level 3 Blas routine.
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
      NOTA = lsame(TRANSA,'N')
      NOTB = lsame(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      ENDIF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      ENDIF
*
*     Test the input parameters.
*
      INFO = 0
      IF ( .NOT.NOTA .AND.  .NOT.lsame(TRANSA,'C') .AND. .
     1      NOT.lsame(TRANSA,'T')) THEN
         INFO = 1
      ELSE IF ( .NOT.NOTB .AND.  .NOT.lsame(TRANSB,'C') .AND.  .NOT.
     1      lsame(TRANSB,'T')) THEN
         INFO = 2
      ELSE IF (M < 0) THEN
         INFO = 3
      ELSE IF (N < 0) THEN
         INFO = 4
      ELSE IF (K < 0) THEN
         INFO = 5
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB < MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC < MAX(1,M)) THEN
         INFO = 13
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO, ' IN GEMM '
         RETURN
      ENDIF
*
*     Quick RETURN IF possible.
*
      IF(M.EQ.0.OR.N.EQ.0.OR.(ALPHA.EQ.ZERO.OR.K.EQ.0).AND.BETA.EQ.ONE)
     1   RETURN
*
*     And IF  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         IF (BETA .EQ. ZERO) THEN
            C(:M,:N) = ZERO
         ELSE
            C(:M,:N) = BETA*C(:M,:N)
         ENDIF
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (NOTB) THEN
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(L,J) .NE. ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(:K,J))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ELSE
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B'' + beta*C
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(J,L) .NE. ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B'' + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(J,:K))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF

      END SUBROUTINE dgemm
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
