c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE slaswp(N, A, LDA, K1, K2, IPIV, INCX)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, LDA, K1, K2, INCX
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(INOUT) :: A
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: I, IP, IX
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL sswap
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, DIMENSION (LDA,N)
!          On ENTRY, the matrix of column DIMENSION N to which the row
!          interchanges will be applied.
!          On EXIT, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, DIMENSION (M*ABS(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF (INCX .EQ. 0) RETURN
      IF (INCX > 0) THEN
         IX = K1
      ELSE
         IX = 1 + (1 - K2)*INCX
      ENDIF
      SELECT CASE (INCX)
      CASE (1)
         DO I = K1, K2
            IP = IPIV(I)
            IF (IP .NE. I) CALL sswap (N, A(I,1), LDA, A(IP,1), LDA)
         END DO
      CASE (2:)
         DO I = K1, K2
            IP = IPIV(IX)
            IF (IP .NE. I) CALL sswap (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      CASE (:(-1))
         DO I = K2, K1, -1
            IP = IPIV(IX)
            IF (IP .NE. I) CALL sswap (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      END SELECT

      END SUBROUTINE slaswp

      SUBROUTINE dlaswp(N, A, LDA, K1, K2, IPIV, INCX)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N, LDA, K1, K2, INCX
      INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
      REAL(WP), DIMENSION(LDA,*), INTENT(INOUT) :: A
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: I, IP, IX
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL dswap
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, DIMENSION (LDA,N)
!          On ENTRY, the matrix of column DIMENSION N to which the row
!          interchanges will be applied.
!          On EXIT, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading DIMENSION of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, DIMENSION (M*ABS(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF (INCX .EQ. 0) RETURN
      IF (INCX > 0) THEN
         IX = K1
      ELSE
         IX = 1 + (1 - K2)*INCX
      ENDIF
      SELECT CASE (INCX)
      CASE (1)
         DO I = K1, K2
            IP = IPIV(I)
            IF (IP .NE. I) CALL dswap (N, A(I,1), LDA, A(IP,1), LDA)
         END DO
      CASE (2:)
         DO I = K1, K2
            IP = IPIV(IX)
            IF (IP .NE. I) CALL dswap (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      CASE (:(-1))
         DO I = K2, K1, -1
            IP = IPIV(IX)
            IF (IP .NE. I) CALL dswap (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      END SELECT

      END SUBROUTINE dlaswp
c !DEC$ ELSE
c       SUBROUTINE la_stub
c       END SUBROUTINE la_stub
c !DEC$ ENDIF

