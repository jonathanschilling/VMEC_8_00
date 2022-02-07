      INTEGER FUNCTION la_ieeeck( ISPEC, ZERO, ONE )
      IMPLICIT NONE
!     GENERIC INTERFACE TO SINGLE OR DOUBLE PRECISION IEEECK ROUTINES
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER            ISPEC
      REAL               ONE, ZERO
C-----------------------------------------------
      INTEGER, EXTERNAL :: ieeeck

      la_ieeeck = ieeeck( ISPEC, ZERO, ONE )

      END FUNCTION la_ieeeck

c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      INTEGER FUNCTION ieeeck( ISPEC, ZERO, ONE )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER            ISPEC
      REAL               ONE, ZERO
C-----------------------------------------------
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic ONLY.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      ieeeck = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         ieeeck = 0
         RETURN
      END IF
*
*
*
*
*     Return IF we were ONLY asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         ieeeck = 0
         RETURN
      END IF

      END FUNCTION ieeeck
c !DEC$ ENDIF
