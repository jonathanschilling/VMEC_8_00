c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE sgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: INFO, LDA, LDB, N, NRHS
      INTEGER :: IPIV(N)
      REAL(WP) :: A(LDA, N), B(LDB, NRHS)
*  purpose
*  =======
*
*  sgesv computes the solution to a real system of linear equations
*     a * x = b,
*  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
*
*  the lu decomposition with partial pivoting and row interchanges is
*  used to factor a as
*     a = p * l * u,
*  where p is a permutation matrix, l is unit lower triangular, and u is
*  upper triangular.  the factored form of a is then used to solve the
*  system of equations a * x = b.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of linear equations, i.e., the order of the
*          matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input/output) real array, dimension (lda,n)
*          on entry, the n-by-n coefficient matrix a.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (output) integer array, dimension (n)
*          the pivot indices that define the permutation matrix p;
*          row i of the matrix was interchanged with row ipiv(i).
*
*  b       (input/output) real array, dimension (ldb,nrhs)
*          on entry, the n-by-nrhs matrix of right hand side matrix b.
*          on exit, if info = 0, the n-by-nrhs solution matrix x.
*
*  ldb     (input) integer
*          the leading dimension of the array b.  ldb >= max(1,n).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero.  the factorization
*                has been completed, but the factor u is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           sgetrf, sgetrs, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL xerbla( 'SGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL sgetrf( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL sgetrs( 'No TRANSPOSE', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF

      END SUBROUTINE sgesv

      SUBROUTINE dgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: INFO, LDA, LDB, N, NRHS
      INTEGER :: IPIV(N)
      REAL(WP) :: A(LDA, N), B(LDB, NRHS)
*  purpose
*  =======
*
*  dgesv computes the solution to a real system of linear equations
*     a * x = b,
*  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
*
*  the lu decomposition with partial pivoting and row interchanges is
*  used to factor a as
*     a = p * l * u,
*  where p is a permutation matrix, l is unit lower triangular, and u is
*  upper triangular.  the factored form of a is then used to solve the
*  system of equations a * x = b.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of linear equations, i.e., the order of the
*          matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input/output) real array, dimension (lda,n)
*          on entry, the n-by-n coefficient matrix a.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (output) integer array, dimension (n)
*          the pivot indices that define the permutation matrix p;
*          row i of the matrix was interchanged with row ipiv(i).
*
*  b       (input/output) real array, dimension (ldb,nrhs)
*          on entry, the n-by-nrhs matrix of right hand side matrix b.
*          on exit, if info = 0, the n-by-nrhs solution matrix x.
*
*  ldb     (input) integer
*          the leading dimension of the array b.  ldb >= max(1,n).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero.  the factorization
*                has been completed, but the factor u is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           dgetrf, dgetrs, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL xerbla( 'dgeSV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL dgetrf( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL dgetrs( 'No TRANSPOSE', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF

      END SUBROUTINE dgesv
c !DEC$ ELSE
c       SUBROUTINE la_stub
c       END SUBROUTINE la_stub
c !DEC$ ENDIF
