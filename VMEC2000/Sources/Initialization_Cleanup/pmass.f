      FUNCTION pmass (xx)
      USE stel_kinds
      USE vmec_input, ONLY: am, bloat,ipmass
      USE vparams, ONLY: mu0
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: xx, pmass, x
      REAL(rprec), parameter :: one=1
C-----------------------------------------------
!     NOTE: On entry, am is in pascals. pmass internal units are mu0*pascals (B**2 units)
      x = MIN (ABS(xx * bloat), 1._dp)

      pmass = 0

      SELECT CASE(ipmass)
      CASE(0)
       DO i = UBOUND(am,1), LBOUND(am,1), -1
        pmass = x*pmass + am(i)
       END DO
      CASE(1)
       pmass = am(0)*(am(1)*(one/(one+(  x/am(2)**2)**am(3))**am(4)
     >                      -one/(one+(one/am(2)**2)**am(3))**am(4))/
     >                  (one-one/(one+(one/am(2)**2)**am(3))**am(4))+
     >          (one-am(1))*(one/(one+(  x/am(5)**2)**am(6))**am(7)
     >                      -one/(one+(one/am(5)**2)**am(6))**am(7))/
     >                  (one-one/(one+(one/am(5)**2)**am(6))**am(7)))
      CASE DEFAULT
        stop 'pmass: Invalid mass profile selected'
      END SELECT
      pmass = mu0*pmass

      END FUNCTION pmass
