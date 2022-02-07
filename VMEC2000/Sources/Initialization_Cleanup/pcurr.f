      FUNCTION pcurr (xx)
      USE stel_kinds
      USE vmec_input, ONLY: ac, bloat, ipcurr
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1        
      INTEGER     :: i, ioff
      REAL(rprec) :: xx, pcurr, x
      REAL(rprec), parameter :: pi=3.141592653589793238462643d0
C-----------------------------------------------
!
!     NOTE:  AC COEFFICIENTS OBTAINED IN THREED1 FILE
!            BY MATCHING TO <JTOR> * dV/dPHI
!
      x = MIN (ABS(xx * bloat), one)

      pcurr = 0

      SELECT CASE(ipcurr)
      CASE(0)
       ioff = LBOUND(ac,1)
       DO i = UBOUND(ac,1), ioff, -1
        pcurr = x*pcurr + ac(i)/(i-ioff+1)
       END DO
       pcurr = x*pcurr
      CASE(1)
       pcurr = ac(1)*(2/pi)*atan(ac(2)*x**ac(3)/(1-x)**ac(4))
     >       + ac(5)*(2/pi)*atan(ac(6)*x**ac(7)/(1-x)**ac(8))
      CASE(2)
       pcurr = ac(1)*(2/pi)*atan(ac(2)*x**ac(3)/(1-x)**ac(4))
     >       + ac(5)*(2/pi)*atan(ac(6)*x**ac(7)/(1-x)**ac(8))
     >       + ac(9)*(2/pi)*atan(ac(10)*x**ac(11)/(1-x)**ac(12))
      CASE DEFAULT
        stop 'pcurr: Invalid mass profile selected'
      END SELECT

      END FUNCTION pcurr
