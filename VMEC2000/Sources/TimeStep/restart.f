      SUBROUTINE restart_iter(time_step)
      USE vmec_main, rprec => rprec1
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec1) :: time_step
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec1), PARAMETER :: c1p03 = 1.03_dp, cp90 = 0.90_dp
c-----------------------------------------------

      SELECT CASE (irst)
      CASE DEFAULT
         xstore(:neqs2) = xc(:neqs2)
         RETURN
      CASE (2:3)
         xcdot(:neqs2) = zero
         xc(:neqs2) = xstore(:neqs2)
         time_step = time_step*((irst-2)/c1p03 + cp90*(3-irst))
         IF (irst .eq. 2) ijacob = ijacob + 1
         irst = 1
         RETURN
      END SELECT

      END SUBROUTINE restart_iter
