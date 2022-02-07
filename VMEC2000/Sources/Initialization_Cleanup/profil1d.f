      SUBROUTINE profil1d(xc, xcdot, lreset)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE vsvd, torflux_edge => torflux
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(neqs2), INTENT(out) :: xc, xcdot
      LOGICAL, INTENT(in) :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: Itor, si, pedge, tflux, tfluxd
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: pcurr, pmass, piota, torflux,
     1    torflux_deriv
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        ai       array of coefficients in phi-series for iota (ncurr=0)
!        ac       array of coefficients in phi-series for the quantity d(Icurv)/ds = toroidal
!                 current density * Vprime, so Icurv(s) = Itor(s) (used for ncurr=1)
!        am       array of coefficients in phi-series for mass (NWT/m**2)
!        iotas    rotational transform , on half radial mesh
!        Icurv    (-)toroidal current inside flux surface (vanishes like s)
!        mass     mass profile on half-grid
!        phiedge  value of REAL toroidal flux at plasma edge (s=1)
!        phips    same as phip , one-dimensional array
!        presf    pressure profile on full-grid, mass/phip**gamma
!        spres_ped value of s beyond which pressure profile is flat (pedestal)

!
!     COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
!     COMPUTE MASS PROFILE ON HALF-GRID
!     BY READING INPUT COEFFICIENTS. PRESSURE CONVERTED TO
!     INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7
!
      torflux_edge = signgs * phifac * phiedge / twopi
      r00 = rmn_bdy(0,0,1)

      phips(1) = 0
      icurv(1) = 0
      icurv(ns+1) = 0         !Need in tomnsps when pushing iota force

      DO i = 2,ns
         si = hs*(i - c1p5)
         tflux = torflux(si)
         phips(i) = torflux_edge * torflux_deriv(si)
         iotas(i) = piota(tflux)
         icurv(i) = pcurr(tflux)
      END DO

      DO i = 1,ns
         si = hs*(i-1)
         tflux = torflux(si)
         iotaf(i) = piota(tflux)
      ENDDO
!
!     SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
!     FACTOR OF SIGNGS NEEDED HERE, SINCE MATCH IS MADE TO LINE
!     INTEGRAL OF BSUBU (IN GETIOTA) ~ SIGNGS * CURTOR
!
      pedge = pcurr(one)
      Itor = signgs*currv/twopi
      IF (ABS(pedge) .gt. ABS(EPSILON(pedge)*curtor))
     1   Itor = Itor/pedge
      icurv(2:ns) = Itor*icurv(2:ns)

!
!     POSSIBLE PRESSURE PEDESTAL FOR S >= SPRES_PED
!
      spres_ped = ABS(spres_ped)
      IF (.not.lrecon) THEN
        DO i = 2,ns
          si = hs*(i - c1p5)
          tflux = torflux(si)
          tfluxd = torflux_edge * torflux_deriv(si)
          IF (si .gt. spres_ped) THEN
             pedge = pmass(spres_ped)
          ELSE
             pedge = pmass(tflux)
          END IF
          mass(i) = pedge*(ABS(tfluxd)*r00)**gamma
        END DO

      ELSE
        iotas(:ns) = 0
        iotaf(:ns) = 0
        mass (:ns) = 0
        presf(:ns) = 0
      END IF

      pres(:ns+1) = 0
      xcdot(:neqs2) = 0

      IF (lreset) THEN
        xc(:neqs1) = 0
        IF (lrecon) iresidue = 0
      END IF
      IF (lrecon) THEN
        IF (iresidue .gt. 1) iresidue = 1
!
!       COMPUTE INDEX ARRAY FOR FINDPHI ROUTINE
!
        DO i = 1,ns
          indexr(i)    = ns + 1 - i                    !FINDPHI
          indexr(i+ns) = i + 1                         !FINDPHI
        ENDDO
        indexr(2*ns) = ns
      END IF

      END SUBROUTINE profil1d
