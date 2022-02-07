      FUNCTION LAMCH_G( CMACH )
      USE LAPREC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      REAL(RPREC) :: LAMCH_G
      CHARACTER(LEN=1), INTENT(IN) :: CMACH
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL, EXTERNAL :: SLAMCH
      DOUBLE PRECISION :: DLAMCH
C-----------------------------------------------
      IF (LDOUBLE) THEN
         LAMCH_G = DLAMCH(CMACH)
      ELSE
         LAMCH_G = SLAMCH(CMACH)
      END IF

      END FUNCTION LAMCH_G

c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      FUNCTION SLAMCH( CMACH )
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      REAL(WP) :: SLAMCH
      CHARACTER(LEN=1), INTENT(IN) :: CMACH
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(WP) :: RMACH
      LOGICAL, EXTERNAL :: la_lsame
*     ..
*
*  Purpose
*  =======
*
*  SLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by SLAMCH:
*          = 'E' or 'e',   SLAMCH := eps
*          = 'S' or 's',   SLAMCH := sfmin
*          = 'B' or 'b',   SLAMCH := base
*          = 'P' or 'p',   SLAMCH := eps*base
*          = 'N' or 'n',   SLAMCH := t
*          = 'R' or 'r',   SLAMCH := rnd
*          = 'M' or 'm',   SLAMCH := emin
*          = 'U' or 'u',   SLAMCH := rmin
*          = 'L' or 'l',   SLAMCH := emax
*          = 'O' or 'o',   SLAMCH := rmax
*
*          WHERE
*
*          eps   = relative machine precision                           EPSILON/RADIX
*          sfmin = safe minimum, such that 1/sfmin does not overflow    TINY
*          base  = base of the machine                                  RADIX
*          prec  = eps*base                                             EPSILON
*          t     = number of (base) digits in the mantissa              DIGITS
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow          MINEXPONENT
*          rmin  = underflow threshold - base**(emin-1)                 TINY
*          emax  = largest exponent before overflow                     MAXEXPONENT
*          rmax  = overflow threshold  - (base**emax)*(1-eps)           HUGE
*
* =====================================================================
*
*     .. Parameters ..
      REAL(WP), PARAMETER :: ONE = 1
*     ..
*     ..
*     .. Executable Statements ..
*
*
      IF( la_lsame( CMACH, 'E' ) ) THEN
         RMACH = EPSILON(ONE)/RADIX(ONE)
      ELSE IF( la_lsame( CMACH, 'S' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( la_lsame( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ONE)
      ELSE IF( la_lsame( CMACH, 'P' ) ) THEN
         RMACH = EPSILON(ONE)
      ELSE IF( la_lsame( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ONE)
      ELSE IF( la_lsame( CMACH, 'R' ) ) THEN
         RMACH = ONE
      ELSE IF( la_lsame( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ONE)
      ELSE IF( la_lsame( CMACH, 'U' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( la_lsame( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ONE)
      ELSE IF( la_lsame( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ONE)
      END IF
*
      SLAMCH = RMACH

      END FUNCTION SLAMCH

      FUNCTION DLAMCH( CMACH )
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      REAL(WP) :: DLAMCH
      CHARACTER(LEN=1), INTENT(IN) :: CMACH
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(WP) :: RMACH
      LOGICAL, EXTERNAL :: la_lsame
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's',   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          WHERE
*
*          eps   = relative machine precision                           EPSILON/RADIX
*          sfmin = safe minimum, such that 1/sfmin does not overflow    TINY
*          base  = base of the machine                                  RADIX
*          prec  = eps*base                                             EPSILON
*          t     = number of (base) digits in the mantissa              DIGITS
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow          MINEXPONENT
*          rmin  = underflow threshold - base**(emin-1)                 TINY
*          emax  = largest exponent before overflow                     MAXEXPONENT
*          rmax  = overflow threshold  - (base**emax)*(1-eps)           HUGE
*
* =====================================================================
*
*     .. Parameters ..
      REAL(WP), PARAMETER :: ONE = 1
*     ..
*     ..
*     .. Executable Statements ..
*
*
      IF( la_lsame( CMACH, 'E' ) ) THEN
         RMACH = EPSILON(ONE)/RADIX(ONE)
      ELSE IF( la_lsame( CMACH, 'S' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( la_lsame( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ONE)
      ELSE IF( la_lsame( CMACH, 'P' ) ) THEN
         RMACH = EPSILON(ONE)
      ELSE IF( la_lsame( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ONE)
      ELSE IF( la_lsame( CMACH, 'R' ) ) THEN
         RMACH = ONE
      ELSE IF( la_lsame( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ONE)
      ELSE IF( la_lsame( CMACH, 'U' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( la_lsame( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ONE)
      ELSE IF( la_lsame( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ONE)
      END IF
*
      DLAMCH = RMACH

      END FUNCTION DLAMCH
c !DEC$ ENDIF
