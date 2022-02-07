c !DEC$ IF DEFINED (NEED_BLAS)
      INTEGER FUNCTION isamax(n, dx, incx)
      USE LAPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      INTEGER, INTENT(IN) :: incx, n
      REAL(WP), INTENT(IN) :: dx(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: i, ix
      REAL(WP) :: dmax
C-----------------------------------------------

      isamax = 0
      IF (n.lt.1 .or. incx.le.0) RETURN
      isamax = 1
      IF (n.eq.1)return
      IF (incx.eq.1) GOTO 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = ABS(dx(1))
      ix = ix + incx
      DO 10 i = 2,n
         if(ABS(dx(ix)).le.dmax) GOTO 5
         isamax = i
         dmax = ABS(dx(ix))
    5    ix = ix + incx
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 dmax = ABS(dx(1))
      DO 30 i = 2,n
         if(ABS(dx(i)).le.dmax) GOTO 30
         isamax = i
         dmax = ABS(dx(i))
   30 CONTINUE

      END FUNCTION isamax

      INTEGER FUNCTION idamax(n, dx, incx)
      USE LAPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      INTEGER, INTENT(IN) :: incx, n
      REAL(WP), INTENT(IN) :: dx(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: i, ix
      REAL(WP) :: dmax
C-----------------------------------------------

      idamax = 0
      IF (n.lt.1 .or. incx.le.0) RETURN
      idamax = 1
      IF (n.eq.1)return
      IF (incx.eq.1) GOTO 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = ABS(dx(1))
      ix = ix + incx
      DO 10 i = 2,n
         if(ABS(dx(ix)).le.dmax) GOTO 5
         idamax = i
         dmax = ABS(dx(ix))
    5    ix = ix + incx
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 dmax = ABS(dx(1))
      DO 30 i = 2,n
         if(ABS(dx(i)).le.dmax) GOTO 30
         idamax = i
         dmax = ABS(dx(i))
   30 CONTINUE

      END FUNCTION idamax
c !DEC$ ELSE
c       SUBROUTINE bla_stub
c       END SUBROUTINE bla_stub
c !DEC$ ENDIF
