c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE sgetrf (M, N, A, LDA, IPIV, INFO)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, LDA
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, DIMENSION(*), INTENT(OUT) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(INOUT) :: A
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: IINFO, J, JB, NB
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL sgemm, sgetf2, slaswp, strsm
      INTEGER, EXTERNAL :: ilaenv
!-----------------------------------------------
!
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  WHERE P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal IF m > n), and U is upper
!  triangular (upper trapezoidal IF m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real array, dimension (LDA,N)
!          On ENTRY, the M-by-N matrix to be factored.
!          On EXIT, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= MAX(1,M).
!
!  IPIV    (output) integer array, dimension (MIN(M,N))
!          The pivot indices; for 1 <= i <= MIN(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
!                has been completed, but the factor u is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M < 0) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (LDA < MAX(1,M)) THEN
         INFO = -4
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN SGETRF'
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     Determine the block size for this environment.
!
      NB = ilaenv(1,'SGETRF',' ',M,N,-1,-1)
!
      IF (NB<=1 .OR. NB>=MIN(M,N)) THEN
!
!        Use unblocked code.
!
         CALL sgetf2 (M, N, A, LDA, IPIV, INFO)
      ELSE
!
!        Use blocked code.
!
         DO J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N) - J + 1,NB)
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL sgetf2 (M - J + 1, JB, A(J,J), LDA, IPIV(J), IINFO)
!
!           Adjust INFO and the pivot indices.
!
            IF (INFO.EQ.0 .AND. IINFO>0) INFO = IINFO + J - 1
            IPIV(J:MIN(M,J+JB-1)) = J - 1 + IPIV(J:MIN(M,J+JB-1))
!
!           Apply interchanges to columns 1:J-1.
!
            CALL slaswp (J - 1, A, LDA, J, J + JB - 1, IPIV, 1)
!
            IF (J + JB <= N) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL slaswp (N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,IPIV,1)
!
!              Compute block row of U. (S/DTRSM) and update trailing
!              submatrix (S/sgemm)
!
               CALL strsm ('Left', 'Lower', 'No TRANSPOSE', 'Unit',
     1                        JB, N - J - JB + 1, ONE, A(J,J),
     2                        LDA, A(J,J+JB), LDA)


               IF (J + JB <= M) CALL sgemm ('No TRANSPOSE',
     1                'No TRANSPOSE', M - J - JB + 1,
     2                 N - J - JB + 1, JB, (-ONE), A(J+JB,J),
     3                 LDA, A(J,J+JB), LDA, ONE, A(J+JB,J+JB), LDA)

            ENDIF
         END DO
      ENDIF

      END SUBROUTINE sgetrf

      SUBROUTINE dgetrf (M, N, A, LDA, IPIV, INFO)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, LDA
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, DIMENSION(*), INTENT(OUT) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(INOUT) :: A
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: IINFO, J, JB, NB
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL dgemm, dgetf2, dlaswp, dtrsm
      INTEGER, EXTERNAL :: ilaenv
!-----------------------------------------------
!
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  dgeTRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  WHERE P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal IF m > n), and U is upper
!  triangular (upper trapezoidal IF m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real array, dimension (LDA,N)
!          On ENTRY, the M-by-N matrix to be factored.
!          On EXIT, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= MAX(1,M).
!
!  IPIV    (output) integer array, dimension (MIN(M,N))
!          The pivot indices; for 1 <= i <= MIN(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
!                has been completed, but the factor u is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M < 0) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (LDA < MAX(1,M)) THEN
         INFO = -4
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN dgeTRF'
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     Determine the block size for this environment.
!
      NB = ilaenv(1,'DGETRF',' ',M,N,-1,-1)
!
      IF (NB<=1 .OR. NB>=MIN(M,N)) THEN
!
!        Use unblocked code.
!
         CALL dgetf2 (M, N, A, LDA, IPIV, INFO)
      ELSE
!
!        Use blocked code.
!
         DO J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N) - J + 1,NB)
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL dgetf2 (M - J + 1, JB, A(J,J), LDA, IPIV(J), IINFO)
!
!           Adjust INFO and the pivot indices.
!
            IF (INFO.EQ.0 .AND. IINFO>0) INFO = IINFO + J - 1
            IPIV(J:MIN(M,J+JB-1)) = J - 1 + IPIV(J:MIN(M,J+JB-1))
!
!           Apply interchanges to columns 1:J-1.
!
            CALL dlaswp (J - 1, A, LDA, J, J + JB - 1, IPIV, 1)
!
            IF (J + JB <= N) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL dlaswp (N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,IPIV,1)
!
!              Compute block row of U. (S/DTRSM) and update trailing
!              submatrix (S/dgemm)
!
               CALL dtrsm ('Left', 'Lower', 'No TRANSPOSE', 'Unit',
     1                        JB, N - J - JB + 1, ONE, A(J,J),
     2                        LDA, A(J,J+JB), LDA)


               IF (J + JB <= M) CALL dgemm ('No TRANSPOSE',
     1                'No TRANSPOSE', M - J - JB + 1,
     2                 N - J - JB + 1, JB, (-ONE), A(J+JB,J),
     3                 LDA, A(J,J+JB), LDA, ONE, A(J+JB,J+JB), LDA)

            ENDIF
         END DO
      ENDIF

      END SUBROUTINE dgetrf
c !DEC$ ELSE
c       SUBROUTINE la_stub
c       END SUBROUTINE la_stub
c !DEC$ ENDIF
