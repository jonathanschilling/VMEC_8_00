      FUNCTION torflux_deriv (x)
      USE stel_kinds
!     USE vmec_input, ONLY: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: x, torflux_deriv
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     torflux_deriv = af(1) + x*(2*af(2) + x*(3*af(3) + x*(4*af(4) +
!    1           x*(5*af(5) + x*(6*af(6) + x*(7*af(7) + x*(8*af(8) +
!    2           x*(9*af(9) + x*10*af(10)))))))))
      torflux_deriv = 1

      END FUNCTION torflux_deriv
