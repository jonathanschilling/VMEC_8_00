      FUNCTION piota (x)
      USE stel_kinds
      USE vmec_input, ONLY: ai, ipiota
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, piota
      REAL(rprec), parameter :: pi=3.141592653589793238462643d0
C-----------------------------------------------
      piota = 0

      SELECT CASE(ipiota)
      CASE(0)
       DO i = UBOUND(ai,1), LBOUND(ai,1), -1
        piota = x*piota + ai(i)
       END DO
      CASE(1)
       piota = ai(0)+ai(1)*(2/pi)*atan(ai(2)*x**ai(3)/(1-x)**ai(4))
     >             + ai(5)*(2/pi)*atan(ai(6)*x**ai(7)/(1-x)**ai(8))
      CASE DEFAULT
       stop 'Invalid mass profile selected'
      END SELECT

      END FUNCTION piota
