c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE sgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => SP
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      INTEGER, INTENT(OUT) :: INFO
      CHARACTER, INTENT(IN) :: TRANS
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB,*), INTENT(INOUT) :: B
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      LOGICAL :: NOTRAN
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL slaswp, strsm
      LOGICAL, EXTERNAL :: lsame
!-----------------------------------------------
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
!  GETRS solves a system of linear equations
!     A * X = B  or  A'' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by GETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No TRANSPOSE)
!          = 'T':  A''* X = B  (Transpose)
!          = 'C':  A''* X = B  (Conjugate TRANSPOSE = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL array, DIMENSION (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by SGETRF.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.  LDA >= MAX(1,N).
!
!  IPIV    (input) INTEGER array, DIMENSION (N)
!          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, DIMENSION (LDB,NRHS)
!          On ENTRY, the right hand side matrix B.
!          On EXIT, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading DIMENSION of the array B.  LDB >= MAX(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful EXIT
!          < 0:  IF INFO = -i, the i-th argument had an illegal value
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
      NOTRAN = lsame(TRANS,'N')
      IF ( .NOT.NOTRAN .AND.  .NOT.lsame(TRANS,'T') .AND.
     1     .NOT.lsame(TRANS,'C')) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (NRHS < 0) THEN
         INFO = -3
      ELSE IF (LDA < MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB < MAX(1,N)) THEN
         INFO = -8
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *, 'INFO = ', INFO, 'IN GETRS'
         RETURN
      ENDIF
!
!     Quick RETURN IF possible
!
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
!

      IF (NOTRAN) THEN
!
!           Solve A * X = B.
!
!           Apply row interchanges to the right hand sides.
            CALL slaswp (NRHS, B, LDB, 1, N, IPIV, 1)
!
!           Solve L*X = B, overwriting B with X.
            CALL strsm('Left', 'Lower', 'No TRANSPOSE', 'Unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
!
!           Solve U*X = B, overwriting B with X.
            CALL strsm('Left', 'Upper', 'No TRANSPOSE', 'Non-unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
      ELSE
!
!           Solve A'' * X = B.
!
!           Solve U''*X = B, overwriting B with X.
            CALL strsm('Left', 'Upper', 'Transpose', 'Non-unit', N,
     1                  NRHS, ONE, A, LDA, B, LDB)
!
!           Solve L''*X = B, overwriting B with X.
            CALL strsm('Left', 'Lower', 'Transpose', 'Unit', N, NRHS,
     1                  ONE, A, LDA, B, LDB)
!
!           Apply row interchanges to the solution vectors.
            CALL slaswp (NRHS, B, LDB, 1, N, IPIV, -1)
      END IF

      END SUBROUTINE sgetrs

      SUBROUTINE dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => DP
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      INTEGER, INTENT(OUT) :: INFO
      CHARACTER, INTENT(IN) :: TRANS
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB,*), INTENT(INOUT) :: B
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      LOGICAL :: NOTRAN
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL dlaswp, dtrsm
      LOGICAL, EXTERNAL :: lsame
!-----------------------------------------------
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
!  GETRS solves a system of linear equations
!     A * X = B  or  A'' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by GETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No TRANSPOSE)
!          = 'T':  A''* X = B  (Transpose)
!          = 'C':  A''* X = B  (Conjugate TRANSPOSE = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL array, DIMENSION (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by SGETRF.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.  LDA >= MAX(1,N).
!
!  IPIV    (input) INTEGER array, DIMENSION (N)
!          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, DIMENSION (LDB,NRHS)
!          On ENTRY, the right hand side matrix B.
!          On EXIT, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading DIMENSION of the array B.  LDB >= MAX(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful EXIT
!          < 0:  IF INFO = -i, the i-th argument had an illegal value
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
      NOTRAN = lsame(TRANS,'N')
      IF ( .NOT.NOTRAN .AND.  .NOT.lsame(TRANS,'T') .AND.
     1     .NOT.lsame(TRANS,'C')) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (NRHS < 0) THEN
         INFO = -3
      ELSE IF (LDA < MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB < MAX(1,N)) THEN
         INFO = -8
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *, 'INFO = ', INFO, 'IN GETRS'
         RETURN
      ENDIF
!
!     Quick RETURN IF possible
!
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
!

      IF (NOTRAN) THEN
!
!           Solve A * X = B.
!
!           Apply row interchanges to the right hand sides.
            CALL dlaswp (NRHS, B, LDB, 1, N, IPIV, 1)
!
!           Solve L*X = B, overwriting B with X.
            CALL dtrsm('Left', 'Lower', 'No TRANSPOSE', 'Unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
!
!           Solve U*X = B, overwriting B with X.
            CALL dtrsm('Left', 'Upper', 'No TRANSPOSE', 'Non-unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
      ELSE
!
!           Solve A'' * X = B.
!
!           Solve U''*X = B, overwriting B with X.
            CALL dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', N,
     1                  NRHS, ONE, A, LDA, B, LDB)
!
!           Solve L''*X = B, overwriting B with X.
            CALL dtrsm('Left', 'Lower', 'Transpose', 'Unit', N, NRHS,
     1                  ONE, A, LDA, B, LDB)
!
!           Apply row interchanges to the solution vectors.
            CALL dlaswp (NRHS, B, LDB, 1, N, IPIV, -1)
      END IF

      END SUBROUTINE dgetrs
c !DEC$ ENDIF

      SUBROUTINE sgetrs1(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => SP
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      INTEGER, INTENT(OUT) :: INFO
      CHARACTER, INTENT(IN) :: TRANS
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB), INTENT(INOUT) :: B

      CALL sgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

      END SUBROUTINE sgetrs1


      SUBROUTINE dgetrs1(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE LAPREC, ONLY: WP => DP
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      INTEGER, INTENT(OUT) :: INFO
      CHARACTER, INTENT(IN) :: TRANS
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(IN) :: A
      REAL(WP), DIMENSION(LDB), INTENT(INOUT) :: B

      CALL dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

      END SUBROUTINE dgetrs1
