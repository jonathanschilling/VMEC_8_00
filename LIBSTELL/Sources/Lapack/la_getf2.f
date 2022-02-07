c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE sgetf2(M, N, A, LDA, IPIV, INFO)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER M, N, LDA, INFO
      INTEGER, DIMENSION(*) :: IPIV
      REAL(WP), DIMENSION(LDA,*) :: A
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: ISMAX (1)
      INTEGER :: J, JP, I
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL sswap, sger
!-----------------------------------------------
!   I N T R I N S I C  F U N C T I O N S
!-----------------------------------------------
      INTRINSIC MAX, MIN
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  WHERE P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal IF m > n), and U is upper
!  triangular (upper trapezoidal IF m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, DIMENSION (LDA,N)
!          On ENTRY, the m by n matrix to be factored.
!          On EXIT, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.  LDA >= MAX(1,M).
!
!  IPIV    (output) INTEGER array, DIMENSION (MIN(M,N))
!          The pivot indices; for 1 <= i <= MIN(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  info    (output) integer
!          = 0: successful exit
!          < 0: if info = -k, the k-th argument had an illegal value
!          > 0: if info = k, u(k,k) is exactly zero. the factorization
!               has been completed, but the factor u is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
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
         PRINT *, 'INFO = ', INFO,' IN SGETF2'
         RETURN
      ENDIF
!
!     Quick RETURN IF possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
      DO J = 1, MIN(M,N)
!
!        find pivot and test for singularity.
!
!        JP = J - 1 + IDAMAX(M - J + 1,A(J,J),1)
         ISMAX = MAXLOC(ABS(A(J:M,J)))
         JP = J - 1 + ISMAX(1)
         IPIV(J) = JP
         IF (A(JP,J) .NE. ZERO) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF (JP .NE. J) CALL sswap (N, A(J,1), LDA, A(JP,1), LDA)
!
!           Compute elements J+1:M of J-th column.
!
!            IF (J < M) CALL DSCAL (M - J, ONE/A(J,J), A(J+1,J), 1)
            A(J+1:M,J) = ONE/A(J,J)*A(J+1:M,J)
!
         ELSE IF (INFO .EQ. 0) THEN
!
            INFO = J
         ENDIF
!
!           IF (J < MIN(M,N)) CALL sger (M - J, N - J, (-ONE), A(J+1,J),
!     1          1, A(J,J+1), LDA, A(J+1,J+1), LDA)


!        Update trailing submatrix.
            IF (J < MIN(M,N)) THEN
                DO i=J+1, N
                   A(J+1:M,i) = A(J+1:M,i) - A(J+1:M,J)*A(J,i)
                END DO
            END IF

      END DO

      END SUBROUTINE sgetf2

      SUBROUTINE dgetf2(M, N, A, LDA, IPIV, INFO)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER M, N, LDA, INFO
      INTEGER, DIMENSION(*) :: IPIV
      REAL(WP), DIMENSION(LDA,*) :: A
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1, ZERO = 0
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: ISMAX (1)
      INTEGER :: J, JP, I
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL dswap, dger
!-----------------------------------------------
!   I N T R I N S I C  F U N C T I O N S
!-----------------------------------------------
      INTRINSIC MAX, MIN
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  dgeTF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  WHERE P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal IF m > n), and U is upper
!  triangular (upper trapezoidal IF m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, DIMENSION (LDA,N)
!          On ENTRY, the m by n matrix to be factored.
!          On EXIT, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.  LDA >= MAX(1,M).
!
!  IPIV    (output) INTEGER array, DIMENSION (MIN(M,N))
!          The pivot indices; for 1 <= i <= MIN(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  info    (output) integer
!          = 0: successful exit
!          < 0: if info = -k, the k-th argument had an illegal value
!          > 0: if info = k, u(k,k) is exactly zero. the factorization
!               has been completed, but the factor u is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
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
         PRINT *, 'INFO = ', INFO,' IN dgeTF2'
         RETURN
      ENDIF
!
!     Quick RETURN IF possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
      DO J = 1, MIN(M,N)
!
!        find pivot and test for singularity.
!
!        JP = J - 1 + IDAMAX(M - J + 1,A(J,J),1)
         ISMAX = MAXLOC(ABS(A(J:M,J)))
         JP = J - 1 + ISMAX(1)
         IPIV(J) = JP
         IF (A(JP,J) .NE. ZERO) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF (JP .NE. J) CALL dswap (N, A(J,1), LDA, A(JP,1), LDA)
!
!           Compute elements J+1:M of J-th column.
!
!            IF (J < M) CALL DSCAL (M - J, ONE/A(J,J), A(J+1,J), 1)
            A(J+1:M,J) = ONE/A(J,J)*A(J+1:M,J)
!
         ELSE IF (INFO .EQ. 0) THEN
!
            INFO = J
         ENDIF
!
!           IF (J < MIN(M,N)) CALL dger (M - J, N - J, (-ONE), A(J+1,J),
!     1          1, A(J,J+1), LDA, A(J+1,J+1), LDA)


!        Update trailing submatrix.
            IF (J < MIN(M,N)) THEN
                DO i=J+1, N
                   A(J+1:M,i) = A(J+1:M,i) - A(J+1:M,J)*A(J,i)
                END DO
            END IF

      END DO

      END SUBROUTINE dgetf2
c !DEC$ ELSE
c       SUBROUTINE la_stub
c       END SUBROUTINE la_stub
c !DEC$ ENDIF
