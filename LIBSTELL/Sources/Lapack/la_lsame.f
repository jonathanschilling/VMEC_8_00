      LOGICAL FUNCTION la_lsame( CA, CB )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      LOGICAL, EXTERNAL :: lsame
C-----------------------------------------------
      la_lsame = lsame( CA, CB )

      END FUNCTION la_lsame

c !DEC$ IF DEFINED (NEED_BLAS)
      LOGICAL FUNCTION lsame( CA, CB )
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. IF CA is the same letter as CB regardless of
*  CASE.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single CHARACTERs to be compared.
*
* =====================================================================
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test IF the CHARACTERs are equal
*
      lsame = CA.EQ.CB
      IF( lsame ) RETURN
*
*     Now test for equivalence IF both CHARACTERs are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper CASE 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper CASE 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     1       INTA.GE.145 .AND. INTA.LE.153 .OR.
     2       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     1       INTB.GE.145 .AND. INTB.LE.153 .OR.
     2       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper CASE 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      lsame = INTA.EQ.INTB
*
      END FUNCTION lsame
c !DEC$ ENDIF
