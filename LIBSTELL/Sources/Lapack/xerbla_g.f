      SUBROUTINE xerbla_g( SRNAME, INFO )
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      CHARACTER*6        SRNAME
      INTEGER            INFO
      CALL xerbla (SRNAME, INFO)

      END SUBROUTINE xerbla_g

c !DEC$ IF DEFINED (NEED_BLAS)
      SUBROUTINE xerbla( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  a message is printed and execution stops.
*
*  installers may consider modifying the stop statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 ) SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', a6, ' parameter number ', i2, ' had ',
     $      'an illegal value' )

      END SUBROUTINE xerbla
c !DEC$ ENDIF
