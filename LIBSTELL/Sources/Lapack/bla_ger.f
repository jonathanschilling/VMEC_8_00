c !DEC$ IF .NOT.DEFINED (SUNOS) .AND. .NOT.DEFINED(OSF1) .AND. .NOT.DEFINED(CRAY)
      SUBROUTINE sger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M,N,INCX,INCY,LDA
      REAL(WP), INTENT(IN) :: ALPHA, X(*), Y(*)
      REAL(WP), INTENT(INOUT) :: A(LDA,*)
!-----------------------------------------------
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y'' + A,
*
*  WHERE alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On ENTRY, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY, ALPHA specifies the scalar alpha.
*           Unchanged on EXIT.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*ABS( INCX ) ).
*           Before ENTRY, the incremented array X must contain the m
*           element vector x.
*           Unchanged on EXIT.
*
*  INCX   - INTEGER.
*           On ENTRY, INCX specifies the increment for the elements of
*           X.
*           Unchanged on EXIT.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*ABS( INCY ) ).
*           Before ENTRY, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on EXIT.
*
*  INCY   - INTEGER.
*           On ENTRY, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before ENTRY, the leading m by n part of the array A must
*           contain the matrix of coefficients. On EXIT, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) program. LDA must be at least MAX(1,m).
*           Unchanged on EXIT.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER :: I,IX,J,JY,KX
      REAL(WP), PARAMETER :: ZERO = 0
      REAL(WP) :: TEMP
      LOGICAL :: OK

      OK = (M.GT.0) .AND. (N.GT.0) .AND. (LDA.GE.M)
*
*
*     Quick RETURN IF possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          DO 20,J = 1,N
             IF (Y(J).NE.ZERO) THEN
                 TEMP = ALPHA*Y(J)
                 DO 10,I = 1,M
                    A(I,J) = A(I,J) + X(I)*TEMP
   10            CONTINUE
             END IF
*
   20     CONTINUE
*
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
*
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN
              JY = 1
*
          ELSE
              JY = 1 - (N-1)*INCY
          END IF
*
          DO 40,J = 1,N
             IF (Y(JY) .NE. ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 IX = KX
                 DO 30,I = 1,M
                    A(I,J) = A(I,J) + X(IX)*TEMP
                    IX = IX + INCX
   30            CONTINUE
             END IF
*
             JY = JY + INCY
   40     CONTINUE
      END IF

      END SUBROUTINE sger

      SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M,N,INCX,INCY,LDA
      REAL(WP), INTENT(IN) :: ALPHA, X(*), Y(*)
      REAL(WP), INTENT(INOUT) :: A(LDA,*)
!-----------------------------------------------
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y'' + A,
*
*  WHERE alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On ENTRY, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on EXIT.
*
*  N      - INTEGER.
*           On ENTRY, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on EXIT.
*
*  ALPHA  - REAL            .
*           On ENTRY, ALPHA specifies the scalar alpha.
*           Unchanged on EXIT.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*ABS( INCX ) ).
*           Before ENTRY, the incremented array X must contain the m
*           element vector x.
*           Unchanged on EXIT.
*
*  INCX   - INTEGER.
*           On ENTRY, INCX specifies the increment for the elements of
*           X.
*           Unchanged on EXIT.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*ABS( INCY ) ).
*           Before ENTRY, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on EXIT.
*
*  INCY   - INTEGER.
*           On ENTRY, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on EXIT.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before ENTRY, the leading m by n part of the array A must
*           contain the matrix of coefficients. On EXIT, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On ENTRY, LDA specifies the first DIMENSION of A as declared
*           in the calling (sub) program. LDA must be at least MAX(1,m).
*           Unchanged on EXIT.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER :: I,IX,J,JY,KX
      REAL(WP), PARAMETER :: ZERO = 0
      REAL(WP) :: TEMP
      LOGICAL :: OK

      OK = (M.GT.0) .AND. (N.GT.0) .AND. (LDA.GE.M)
*
*
*     Quick RETURN IF possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          DO 20,J = 1,N
             IF (Y(J).NE.ZERO) THEN
                 TEMP = ALPHA*Y(J)
                 DO 10,I = 1,M
                    A(I,J) = A(I,J) + X(I)*TEMP
   10            CONTINUE
             END IF
*
   20     CONTINUE
*
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
*
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN
              JY = 1
*
          ELSE
              JY = 1 - (N-1)*INCY
          END IF
*
          DO 40,J = 1,N
             IF (Y(JY) .NE. ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 IX = KX
                 DO 30,I = 1,M
                    A(I,J) = A(I,J) + X(IX)*TEMP
                    IX = IX + INCX
   30            CONTINUE
             END IF
*
             JY = JY + INCY
   40     CONTINUE
      END IF

      END SUBROUTINE dger
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
