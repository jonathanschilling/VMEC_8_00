c !DEC$ IF DEFINED (NEED_BLAS)
      FUNCTION snrm2 (N, SX, INCX)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER N, INCX
      REAL(WP), DIMENSION(N) :: SX
      REAL(WP) :: snrm2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
!     REAL(WP), PARAMETER :: CUTLO=4.441E-16_DP, CUTHI=1.304E19_DP
      INTEGER :: NEXT, NN, I, J
      REAL(WP) :: CUTLO, CUTHI, HITEST, NSUM, XMAX
C-----------------------------------------------
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C
      CUTLO = SQRT(TINY(CUTLO)/EPSILON(CUTLO))
      CUTHI = SQRT(HUGE(CUTHI))

      IF (N <= 0) THEN
         snrm2 = ZERO
      ELSE
C
         NEXT = 30
         NSUM = ZERO
         NN = N*INCX
C                                                 BEGIN MAIN LOOP
         I = 1
   20    CONTINUE
         IF (NEXT == 30) GO TO 30
         IF (NEXT == 50) GO TO 50
         IF (NEXT == 70) GO TO 70
         IF (NEXT == 110) GO TO 110
   30    CONTINUE
         IF (ABS(SX(I)) > CUTLO) GO TO 85
         NEXT = 50
         XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50    CONTINUE
         IF (SX(I) == ZERO) GO TO 200
         IF (ABS(SX(I)) > CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
         NEXT = 70
         GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100    CONTINUE
         I = J
         NEXT = 110
         NSUM = (NSUM/SX(I))/SX(I)
  105    CONTINUE
         XMAX = ABS(SX(I))
         GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70    CONTINUE
         IF (ABS(SX(I)) > CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110    CONTINUE
         IF (ABS(SX(I)) <= XMAX) GO TO 115
         NSUM = ONE + NSUM*(XMAX/SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
C
  115    CONTINUE
         NSUM = NSUM + (SX(I)/XMAX)**2
         GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75    CONTINUE
         NSUM = (NSUM*XMAX)*XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85    CONTINUE
         HITEST = CUTHI/N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
         DO J = I, NN, INCX
            IF (ABS(SX(J)) >= HITEST) GO TO 100
            NSUM = NSUM + SX(J)**2
         END DO
         snrm2 = SQRT(NSUM)
         GO TO 300
C
  200    CONTINUE
         I = I + INCX
         IF (I <= NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
         snrm2 = XMAX*SQRT(NSUM)
      ENDIF
  300 CONTINUE

      END FUNCTION snrm2

      FUNCTION dnrm2 (N, SX, INCX)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER N, INCX
      REAL(WP), DIMENSION(N) :: SX
      REAL(WP) :: dnrm2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
!     REAL(WP), PARAMETER :: CUTLO=4.441E-16_DP, CUTHI=1.304E19_DP
      INTEGER :: NEXT, NN, I, J
      REAL(WP) :: CUTLO, CUTHI, HITEST, NSUM, XMAX
C-----------------------------------------------
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C
      CUTLO = SQRT(TINY(CUTLO)/EPSILON(CUTLO))
      CUTHI = SQRT(HUGE(CUTHI))

      IF (N <= 0) THEN
         dnrm2 = ZERO
      ELSE
C
         NEXT = 30
         NSUM = ZERO
         NN = N*INCX
C                                                 BEGIN MAIN LOOP
         I = 1
   20    CONTINUE
         IF (NEXT == 30) GO TO 30
         IF (NEXT == 50) GO TO 50
         IF (NEXT == 70) GO TO 70
         IF (NEXT == 110) GO TO 110
   30    CONTINUE
         IF (ABS(SX(I)) > CUTLO) GO TO 85
         NEXT = 50
         XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50    CONTINUE
         IF (SX(I) == ZERO) GO TO 200
         IF (ABS(SX(I)) > CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
         NEXT = 70
         GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100    CONTINUE
         I = J
         NEXT = 110
         NSUM = (NSUM/SX(I))/SX(I)
  105    CONTINUE
         XMAX = ABS(SX(I))
         GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70    CONTINUE
         IF (ABS(SX(I)) > CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110    CONTINUE
         IF (ABS(SX(I)) <= XMAX) GO TO 115
         NSUM = ONE + NSUM*(XMAX/SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
C
  115    CONTINUE
         NSUM = NSUM + (SX(I)/XMAX)**2
         GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75    CONTINUE
         NSUM = (NSUM*XMAX)*XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85    CONTINUE
         HITEST = CUTHI/N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
         DO J = I, NN, INCX
            IF (ABS(SX(J)) >= HITEST) GO TO 100
            NSUM = NSUM + SX(J)**2
         END DO
         dnrm2 = SQRT(NSUM)
         GO TO 300
C
  200    CONTINUE
         I = I + INCX
         IF (I <= NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
         dnrm2 = XMAX*SQRT(NSUM)
      ENDIF
  300 CONTINUE

      END FUNCTION dnrm2
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF

