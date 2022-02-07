      FUNCTION torflux (x)
      USE stel_kinds
!     USE vmec_input, ONLY: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: x, torflux
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     NOTE: af are normed so that SUM(af) = 1, so that torflux(1) = 1
!
!     torflux = x*(af(1) + x*(af(2) + x*(af(3) + x*(af(4) +
!    1     x*(af(5) + x*(af(6) + x*(af(7) + x*(af(8) + x*(af(9) +
!    2     x*af(10))))))))))
      torflux = x

      END FUNCTION torflux
