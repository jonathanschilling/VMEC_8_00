      SUBROUTINE eqfor(br, bz, bsubu, bsubv, tau, rzl_array)
      USE vmec_main
      USE vmec_params
      USE realspace
      USE vforces, r12 => armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,
     2   bphi => czmn_o
      USE vacmod
      USE vsvd
      USE vspline
      USE csplinx
      USE vmec_io
      USE mgrid_mod
      USE fbal
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt,0:1), INTENT(in) :: bsubu, bsubv
      REAL(rprec), DIMENSION(nrzt), INTENT(out) :: br, bz
      REAL(rprec), DIMENSION(nrzt), INTENT(out) :: tau
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1  TARGET, INTENT(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, icount, itheta, js, js1, l, loff,
     1   lpi, lt, n, n1, nchicur, nchiiota0, ncol, noff,
     2   nout, nsort, iv, iu, lk, nplanes
      REAL(rprec), DIMENSION(:), POINTER ::
     1   rmags, zmags, rmaga, zmaga
      REAL(rprec), DIMENSION(:,:,:), POINTER :: rmncc,zmnsc
      REAL(rprec), DIMENSION(ns) :: phi1, chi1, bucof, bvcof, jPS2
      REAL(rprec) :: modb(nznt)
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   bsup1u, dbts1u, dint1u, t12u, guu_1u, guus1u, r3v,
     2   redg1u, rbps1u
      REAL(rprec) :: aminr1, aminr2, aminr2in, anorm,
     1   aspectratio, betai, betstr,
     2   bminz2, bminz2in, btor, iotamax, musubi,
     3   bzcalc, bzin, chisq, chiwgt, cur0,
     4   delphid_exact, delta1, delta2, delta3, denwgt, lambda,
     5   dlogidsi, dmusubi_meas, er, es, fac, facnorm, factor, fgeo,
     6   fmax, fmin, flao, fpsi0, pavg, pitchc, pitchm,
     7   pprime, qedge, qmin1, qmin2, qmin3, qzero,
     8   raxis0, rcalc, rcen, rcenin, rgeo, rjs, 
     9   rjs1, rlao, rqmin1, rqmin2, rshaf, rshaf1, rshaf2, s11, s12,
     A   s13, s2, s3, sigr0, sigr1, sigz1, smaleli,
     B   splintx, splints, sqmin, sumbpol, sumbtot, sumbtor, sump,
     C   sump2, sump20, t1, tz, jpar_perp=0, jparPs_perp=0,
     D   tol, toroidal_flux, vnorm, vprime, wght0, xmax,
     E   xmida, xmidb, xmin, rzmax, rzmin, zxmax, zxmin, zaxis0,
     F   zmax, zmin, yr1u, yz1u, waist(2), height(2)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
       EXTERNAL splintx,splints
C-----------------------------------------------
!
!     POINTER ASSOCIATIONS
!
      rmags => rzl_array(1,:,0,rcc)
      zmags => rzl_array(1,:,0,zcs+ntmax)
      rmncc => rzl_array(:,:,:,rcc)
      zmnsc => rzl_array(:,:,:,zsc+ntmax) 
      IF (lasym) THEN
        rmaga => rzl_array(1,:,0,rsc)
        zmaga => rzl_array(1,:,0,zcc+ntmax)
      END IF


      CALL bss (r12, bzmn, brmn, azmn, armn, crmn_o, bsupu, bsupv,
     1          br, bphi, bz)

!
!     STORE EDGE VALUES OF B-FIELD
!
      IF (lfreeb .or. ledge_dump) THEN
         IF (ALLOCATED(bredge)) DEALLOCATE (bredge, bpedge, bzedge)
         ALLOCATE (bredge(2*nznt), bpedge(2*nznt), bzedge(2*nznt),
     1             stat=i)
         IF (i .ne. 0) STOP 'Error in EQFOR allocatin bredge'
         DO iv = 1,nzeta
            DO iu = 1,ntheta3
               lk = iv + nzeta*(iu-1)
               n1 = ns*lk
               bredge(lk) = 1.5_dp*br(n1)   - p5*br(n1-1)
               bpedge(lk) = 1.5_dp*bphi(n1)   - p5*bphi(n1-1)
               bzedge(lk) = 1.5_dp*bz(n1)   - p5*bz(n1-1)
            END DO
         END DO
      END IF

!
!
!     NOTE: JXBFORCE ROUTINE MUST BE CALLED TO COMPUTE IZETA, JDOTB
!           ON OUTPUT, J, IZETA, JDOTB ARE IN MKS UNITS (HAVE 1/MU0)
!
!     CAUTION: THIS CALL WILL WRITE OVER br, bz 
!
      CALL jxbforce (bsupu, bsupv, bsubu, bsubv, crmn_o,
     1               rcon, zcon, rcon0, zcon0, gsqrt, izeta, bsq)

!
!     HALF-MESH VOLUME-AVERAGED BETA
!
      tau(1) = 0
      tau(2:nrzt) = signgs*wint(2:nrzt)*gsqrt(2:nrzt)
      DO i = 2, ns
         s2 = SUM(bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i) - pres(i)
         overr(i) = SUM(tau(i:nrzt:ns)/r12(i:nrzt:ns)) / vp(i)
         beta_vol(i) = pres(i)/s2
      END DO

      betaxis = c1p5*beta_vol(2) - p5*beta_vol(3)


      WRITE (nthreed, 5)
    5 FORMAT(/,' NOTE: <RADIAL FORCE> = d(Ipol)/dPHI',
     1         ' - IOTA*d(Itor)/dPHI - dp/dPHI * dV/dPHI',/,
     1         '                      = dV/dPHI*[<JSUPU>',
     1         ' - IOTA*<JSUPV> - SIGN(JAC)*dp/dPHI]',/,
     1         ' (NORMED TO SUM OF INDIVIDUAL TERMS)',//,
     1         '      S     <RADIAL   TOROIDAL     IOTA     ',
     1         ' <JSUPU>    <JSUPV>     d(VOL)/',
     2         '   d(PRES)/    <M>     PRESF    <BSUBU>    <BSUBV>',
     3         '      <J.B>      <B.B>',/,
     4         '             FORCE>      FLUX                ',
     5         '                       d(PHI) ',
     6         '    d(PHI)                             ',/,148('-'),/)

      phipf(1) = twopi*signgs*(c1p5*phip(2) - p5*phip(3))
      presf(1) = c1p5*pres(2) - p5*pres(3)
      DO i = 2,ns1
         presf(i) = p5*(pres(i) + pres(i+1))
         phipf(i) = p5*twopi*signgs*(phip(i) + phip(i+1))
      END DO
      presf(ns) = c1p5*pres(ns)- p5*pres(ns-1)
      phipf(ns) = twopi*signgs*(c1p5*phip(ns) - p5*phip(ns1))

      phi1(1) = zero
      chi1(1) = zero
      DO i = 2, ns
         phi1(i) = phi1(i-1) + hs*phip(i)
         chi1(i) = chi1(i-1) + hs*(phip(i)*iotas(i))
      END DO

      CALL calc_fbal(bsubu, bsubv)

      bucof(1) = 0
      bvcof(1) = c1p5*bvco(2) - p5*bvco(3)
!
!     NOTE:  jcuru, jcurv on FULL radial mesh coming out of calc_fbal
!            They are local (surface-averaged) current densities (NOT integrated in s)
!
      DO i = 2,ns1
         equif(i) = equif(i)*vpphi(i)/(ABS(jcurv(i)*iotaf(i))
     1            + ABS(jcuru(i))+ABS((presgrad(i)*vpphi(i))))
         bucof(i) = p5*(buco(i) + buco(i+1))
         bvcof(i) = p5*(bvco(i) + bvco(i+1))
      END DO

      bucof(ns) = c1p5*buco(ns) - p5*buco(ns1)
      bvcof(ns) = c1p5*bvco(ns) - p5*bvco(ns1)

      equif(1) = two*equif(2) - equif(3)
      jcuru(1) = two*jcuru(2) - jcuru(3)
      jcurv(1) = two*jcurv(2) - jcurv(3)
      presgrad(1)  = two*presgrad(2) - presgrad(3)
      presgrad(ns) = two*presgrad(ns1) - presgrad(ns1-1)
      vpphi(1)  = two*vpphi(2) - vpphi(3)
      vpphi(ns) = two*vpphi(ns1) - vpphi(ns1-1)
      equif(ns) = two*equif(ns1) - equif(ns1-1)
      jcuru(ns) = two*jcuru(ns1) - jcuru(ns1-1)
      jcurv(ns) = two*jcurv(ns1) - jcurv(ns1-1)
      fac = twopi*signgs
      DO js = 1, ns
         es = (js - 1)*hs
         cur0 = signgs*fac*vpphi(js)*phipf(js)
         WRITE (nthreed, 30) es, equif(js), fac*phi1(js),
     1     iotaf(js), fac*jcuru(js)/cur0/mu0, fac*jcurv(js)/cur0/mu0,
     2     fac*vpphi(js), presgrad(js)/phipf(js)/mu0, specw(js), 
     3     presf(js)/mu0, bucof(js), bvcof(js), jdotb(js), bdotb(js)
      END DO
 30   FORMAT(1p,2e10.2,6e11.3,0p,f7.3,1p,5e11.3)

!
!     Calculate mean (toroidally averaged) poloidal cross section area & toroidal flux
!
      anorm = twopi*hs
      vnorm = twopi*anorm
      toroidal_flux = anorm * SUM(bsupv(2:nrzt)*tau(2:nrzt))

!
!     Calculate poloidal circumference and normal surface area and aspect ratio
!     Normal is | dr/du X dr/dv | = SQRT [R**2 guu + (RuZv - RvZu)**2]
!
      ALLOCATE (guu_1u(nznt), guus1u(nznt))
      guu_1u(:nznt) = ru0(ns:nrzt:ns)*ru0(ns:nrzt:ns) +
     1   zu0(ns:nrzt:ns)*zu0(ns:nrzt:ns)
      guus1u(:nznt) = wint(ns:nrzt:ns)*SQRT(guu_1u(:nznt))
      circum_p = twopi*SUM(guus1u(:nznt))
      guus1u(:nznt) = wint(ns:nrzt:ns)*SQRT(
     1     + (r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1))**2*guu_1u(:nznt)
     2     +((rv(ns:nrzt:ns,0) + rv(ns:nrzt:ns,1))*zu0(ns:nrzt:ns) 
     3     - (zv(ns:nrzt:ns,0) + zv(ns:nrzt:ns,1))*ru0(ns:nrzt:ns))**2)
      surf_area_p = twopi**2*SUM(guus1u(:nznt))
      DEALLOCATE (guu_1u, guus1u)

      aspect = aspectratio()

!     Also, estimate mean elongation of plasma from the following relations
!     
!     surf_area_p = 2*pi*R * 2*pi*a sqrt(1+kappa_p**2)
!     volume_p    = 2*pi*R * pi*a**2 * kappa_p
!     cross_area_p=   pi*a**2 * kappa_p
!

      aminr1 = 2*volume_p/surf_area_p
      er = aminr1/aminor_p
      er = one/er**2
      kappa_p = er + SQRT(ABS(er*er - 1))
 
!
!
!     OUTPUT BETAS, INDUCTANCES, SAFETY FACTORS, ETC.
!     (EXTRACTED FROM FQ-CODE, 9-10-92)
!
!     b poloidals (cylindrical estimates)
!
      rcen = p5*(router + rinner)               !geometric center
      n = 0
      n1 = n + 1
      rcenin = DOT_PRODUCT(rmncc(ns,n1,:mpol1+1:2),
     1                     mscale(:mpol1:2)*nscale(n))

      l = (mpol1+1)/2
      ALLOCATE (t12u(l))
      t12u(:l) = mscale(1:mpol1:2)*nscale(n)
      aminr2in = DOT_PRODUCT(rmncc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2in = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2 = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      DEALLOCATE (t12u)
      aminr1 = SQRT(two*volume_p/(twopi*twopi*r00))  !vol av minor radius
      aminr2 = p5*(router - rinner)             !geometric plasma radius
!
!       cylindrical estimates for beta poloidal
      sump = vnorm*SUM(vp(2:ns)*pres(2:ns))
      pavg = sump/volume_p
      ppeak = presf(1)
      factor = two*pavg
!
!       delphid_exact = Integral[ (Bvac - B) * dSphi ]
!
      delphid_exact = zero                         !Eq. 20 in Shafranov
      musubi = zero

      rshaf1 = zero
      rshaf2 = zero
      ALLOCATE (bsup1u(nznt), dbts1u(nznt), dint1u(nznt))
      DO js = 2, ns
         bsup1u(:nznt) = bsubvvac/(r12(js:nrzt:ns)*
     1      r12(js:nrzt:ns))
         delphid_exact = delphid_exact + SUM((bsup1u(:nznt) -
     1      bsupv(js:nrzt:ns))*tau(js:nrzt:ns))
         dbts1u(:nznt) = bsup1u(:nznt)*bsubvvac -
     1      bsupv(js:nrzt:ns)*bsubv(js,:nznt,0)
         musubi = musubi + SUM(dbts1u(:nznt)*tau(js:nrzt:ns))
         dint1u(:nznt) = (bsupu(js:nrzt:ns)*bsubu(js,:nznt,0)
     1      + dbts1u(:nznt) + two*pres(js))*tau(js:nrzt:ns)
         rshaf1 = rshaf1 + SUM(dint1u(:nznt))
         rshaf2 = rshaf2 + SUM(dint1u(:nznt)/r12(js:nrzt:ns))
      END DO
      DEALLOCATE (bsup1u, dbts1u, dint1u)

      delphid_exact = anorm*delphid_exact
      rshaf = rshaf1/rshaf2
      fpsi0 = c1p5*bvco(2) - p5*bvco(3)
      b0 = fpsi0/r00

      IF (lrecon) THEN
      IF (apres .ne. aminr1) THEN
         WRITE (*, 50) (apres/aminr1)**2
         WRITE (nthreed, 50) (apres/aminr1)**2
      ENDIF
   50 FORMAT(/'  Multiply final PHIEDGE by ',f7.3,
     1   ' to make apres = aminr1')
!
!     Output some of measured data
!

      raxis0 = SUM(raxis_cc(0:ntor))
      zaxis0 = SUM(zaxis_cs(0:ntor))
      cur0 = signgs*iotaf(1)*fpsi0/r00**2*(dkappa + one/dkappa)
      IF (lpprof) WRITE (nthreed, 60) presfac*pfac, (-100.*(r00 -
     1   rthompeak))
      WRITE (nthreed, 65) b0, cur0/mu0
   60 FORMAT(//,' Input pressure scaled by ',f6.3/,
     1   ' Pressure profile shifted ',f6.2,' cm.',
     2   ' relative to magnetic axis')
   65 FORMAT(//' B-PHI(R=Raxis,Z=0) = ',f8.3,' [Wb/M**2]',/,
     1   ' J-PHI(R=Raxis,Z=0) = iota(0)*Bt0/(mu0*R0)*(1+k**2)/k = ',1p
     2   e10.3,' [A/M**2]',2/,' Comparison of Input vs Output',8x,
     3   'Input',9x,'VMEC Output',/,1x,71('-'))
      WRITE (nthreed, 70) 1.e-6_dp*currv/mu0, 1.e-6_dp*ctor/mu0,
     1   phiedge, toroidal_flux, 
     1   1.e3_dp*phidiam, 1.e3_dp*delphid, 1.e-3_dp*pthommax,
     2   1.e-3_dp*ppeak/mu0, rcenin, rcen, raxis0, r00, zaxis0, z00,
     2   rc0mse,
     3   apres, aminr1, aminr2in, aminr2, bminz2in, bminz2, rinner,
     4   router, 1.e3_dp*delphid_exact
   70 FORMAT(' Toroidal Current     =',2(10x,f10.3),'  [MA]',/,
     1   ' Edge Toroidal Flux   =',2(10x,f10.3),'  [Wb]',/,
     2   ' Diamagnetic Flux(*)  =',2(10x,f10.3),'  [mWb]',/,
     3   ' Peak Pressure        =',2(10x,f10.3),'  [KPa]',/,
     4   ' Geometric Center     =',2(10x,f10.3),'  [M]',/,
     5   ' Magnetic R Axis      =',2(10x,f10.3),'  [M]',/,
     6   ' Magnetic Z Axis      =',2(10x,f10.3),'  [M]',/,
     7   ' MSE R Axis           =',10x,f10.3,20x,'  [M]',/,
     8   ' Minor Radius (apres) =',2(10x,f10.3),'  [M]',/,
     9   ' Minor Radius (a)     =',2(10x,f10.3),'  [M]',/,
     .   ' Minor Radius (b)     =',2(10x,f10.3),'  [M]',/,
     1   ' Inboard  Midplane R  =',30x,f10.3,'  [M]',/,
     2   ' Outboard Midplane R  =',30x,f10.3,'  [M]',/,50('-')/,
     3   ' * Exact diamagnetic flux (not based on equilibrium) = ',f8.3,
     4   '  [mWb]')

!       Calculate components of total chi-squared
      total_chi_cur = (currv - ctor)**2/sigma_current**2
      nchicur = 1
      IF (iphidiam .eq. 1) total_chi_dia = (phidiam - delphid)**2/
     1   sigma_delphid**2
      nchidia = 1
      END IF   !!IF(lrecon)

      rmax_surf = MAXVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      rmin_surf = MINVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      zmax_surf = MAXVAL(ABS(z1(ns:nrzt:ns,0)+z1(ns:nrzt:ns,1)))

      DO js = 2, ns
         modb(:nznt) = SQRT(two*(bsq(js:nrzt:ns)-pres(js)))
         CALL bextrema (modb, bmin(1,js), bmax(1,js), nzeta, ntheta2)
      END DO

!
!     output geometrical, |B| quantities
!
      CALL elongation (r1, z1, waist, height)

      WRITE (nthreed, 75) bmin(1,ns), bmax(1,ns), bmin(ntheta2,ns), bmax
     1   (ntheta2,ns)
   75 FORMAT(/
     1   ' Magnetic field modulation (averaged over toroidal angle)',/,
     2   1x,71('-')/,' BMIN(u=0)             = ',f14.6/
     3   ' BMAX(u=0)             = ',f14.6/' BMIN(u=pi)            = ',
     4   f14.6/' BMAX(u=pi)            = ',f14.6/)

      sumbtot = 2*(vnorm*SUM(bsq(2:nrzt)*tau(2:nrzt)) - sump)
      sumbtor = vnorm*SUM(tau(2:nrzt)*(r12(2:nrzt)*bsupv(2:nrzt))**2) 
      sumbpol = sumbtot - sumbtor
      betapol = two*sump/sumbpol
      sump20 = two*sump
      sump2 = zero
      sump2 = SUM(pres(2:ns)*pres(2:ns)*vp(2:ns)*vnorm)
      betatot = sump20/sumbtot
      betator = sump20/sumbtor
      VolAvgB = SQRT(ABS(sumbtot/volume_p))
      IonLarmor = 0.0032_dp/VolAvgB
      jPS2(2:ns1) = (jpar2(2:ns1) - jdotb(2:ns1)**2/bdotb(2:ns1))
      jpar_perp = SUM(jpar2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      jparPS_perp = SUM(jPS2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      s2 = SUM(jperp2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      IF (s2 .ne. zero) THEN
         jpar_perp = jpar_perp/s2
         jparPS_perp = jparPS_perp/s2
      END IF
      IF (ntor .gt. 1) THEN
      WRITE (nthreed, 80) aspect, kappa_p, volume_p, cross_area_p, 
     1   surf_area_p, circum_p, Rmajor_p, Aminor_p, rmin_surf, 
     2   rmax_surf, zmax_surf, waist(1), height(1), waist(2), height(2)
      ELSE
      WRITE (nthreed, 80) aspect, kappa_p, volume_p, cross_area_p, 
     1   surf_area_p, circum_p, Rmajor_p, Aminor_p, rmin_surf, 
     2   rmax_surf, zmax_surf, waist(1), height(1)
      END IF
 80   FORMAT(/,' Geometric and Magnetic Quantities',/,1x,71('-')/,
     1   ' Aspect Ratio          = ',f14.6, /
     1   ' Mean Elongation       = ',f14.6, /
     1   ' Plasma Volume         = ',f14.6,' [M**3]',/
     2   ' Cross Sectional Area  = ',f14.6,' [M**2]',/
     3   ' Normal Surface Area   = ',f14.6,' [M**2]',/
     4   ' Poloidal Circumference= ',f14.6,' [M]',/
     5   ' Major Radius          = ',f14.6,' [M]',
     6   ' (from Volume and Cross Section)',/
     7   ' Minor Radius          = ',f14.6,' [M]',
     8   ' (from Cross Section)',/
     9   ' Minimum (inboard)  R  = ',f14.6,' [M]',/
     A   ' Maximum (outboard) R  = ',f14.6,' [M]',/
     A   ' Maximum height     Z  = ',f14.6,' [M]',/
     B   ' Waist (v = 0)   in R  = ',f14.6,' [M]',/
     B   ' Full Height(v = 0)    = ',f14.6,' [M]',:,/
     B   ' Waist (v = pi)  in R  = ',f14.6,' [M]',:,/
     B   ' Full Height(v = pi)   = ',f14.6,' [M]')
      WRITE (nthreed, 85) toroidal_flux, 1.e-6_dp*ctor/mu0, rbtor, 
     1       rbtor0, VolAvgB, IonLarmor, jpar_perp, jparPS_perp
 85   FORMAT(
     1   ' Toroidal Flux         = ',f14.6,' [Wb]',/
     1   ' Toroidal Current      = ',f14.6,' [MA]',/
     1   ' RBtor(s=1)            = ',f14.6,' [T-m]',/
     1   ' RBtor(s=0)            = ',f14.6,' [T-m]',/
     1   ' Volume Average B      = ',f14.6,' [T]',/
     2   ' Ion Larmor Radius     = ',f14.6,' [M] X Ti(keV)**0.5',/
     3   ' <J||**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/
     4   ' <JPS**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/)

      WRITE (nthreed, 90)
   90 FORMAT(/,71('-'),/,' MORE GEOMETRIC AND PHYSICS QUANTITIES',/,
     1    71('-'),/,' Toroidal Plane: Phi = 0',/,
     1    5x,'j',3x,'psi-psiaxis',9x,'a [M]',3x,'ellipticity',3x,
     2   'indentation',7x,'d-shape',4x,'rel. shift',6x,'<J||**2>/',4x,
     3   '<JPS**2>/',/,95x,
     4   '<J-perp**2>',3x,'<J-perp**2>'/,' -----',8(2x,12('-')))

      fac = twopi*hs*signgs
      psi(1) = zero
      ALLOCATE (r3v(ns-1))
      r3v(:ns-1) = fac*phip(2:ns)*iotas(2:ns)
      DO i = 1, ns - 1
         psi(1+i) = psi(i) + r3v(i)
      END DO
      DEALLOCATE (r3v)

!     nphi-plane, noff = 1,....,nzeta
      DO nplanes = 1, 2
         IF (nplanes .eq. 1) THEN
            noff = 1
         ELSE
            IF (nzeta .eq. 1) EXIT
            WRITE (nthreed, 95)
            noff = 1 + nzeta/2
         END IF
          
         ygeo(1) = zero
         DO js = 2, ns
            zmin =  HUGE(zmin)
            zmax = -HUGE(zmax)
            xmin =  HUGE(xmin)
            xmax = -HUGE(xmax)
            rzmax = zero

!           Theta = 0 to pi in upper half of X-Z plane
            DO icount = 1,2
               n1 = noff                                      !nphi-plane, n1 = noff,...,nzeta
               IF (icount .eq. 2)
     1         n1 = MOD(nzeta + 1 - noff,nzeta) + 1           !(twopi-v), reflected plane
               loff = js + ns*(n1-1)
               t1 = one
               IF (icount .eq. 2) t1 = -one
               DO itheta = 1,ntheta2
                  yr1u = r1(loff,0) + sqrts(js)*r1(loff,1)
                  yz1u = z1(loff,0) + sqrts(js)*z1(loff,1)
                  yz1u = t1*yz1u
                  IF (yz1u .ge. zmax) THEN
                     zmax = ABS(yz1u)
                     rzmax = yr1u
                  ELSE IF (yz1u .le. zmin) THEN
                     zmin = yz1u
                     rzmin = yr1u
                  END IF
                  IF (yr1u .ge. xmax) THEN
                     xmax = yr1u
                     zxmax = yz1u
                  ELSE IF (yr1u .le. xmin) THEN
                     xmin = yr1u
                     zxmin = yz1u
                  END IF
                  loff = loff + ns*nzeta
               END DO
            END DO


            lpi = ns*((noff-1) + nzeta*(ntheta2 - 1))
            lt  = ns*((noff-1))
            xmida = r1(js+lpi,0) + sqrts(js)*r1(js+lpi,1)
            xmidb = r1(js+lt,0)  + sqrts(js)*r1(js+lt,1)

            rgeo = p5*(xmidb + xmida)              !Geometric major radius
            ygeo(js) = p5*(xmidb - xmida)          !Geometric minor radius
   
            yinden(js) = (xmida - xmin)/(xmax - xmin) !Geometric indentation
            yellip(js) = (zmax - zmin)/(xmax - xmin)  !Geometric ellipticity
   
            ytrian(js) = (rgeo - rzmax)/(xmax - xmin) !Geometric triangularity
            yshift(js) = (r1(1+lt,0)-rgeo)/(xmax - xmin) !Geometric shift

            IF (jperp2(js) .eq. zero) jperp2(js) = EPSILON(jperp2(js))
            jpar_perp = jpar2(js)/jperp2(js)
            IF (js .lt. ns) THEN
               jparPS_perp = jPS2(js)/jperp2(js)
            ELSE
               jparPS_perp = zero
            END IF

            IF (nplanes .eq. 1) THEN
               WRITE (nthreed, 120) js, psi(js), ygeo(js), yellip(js),
     1            yinden(js), ytrian(js), yshift(js), jpar_perp, 
     2            jparPS_perp
            ELSE
               WRITE (nthreed, 120) js, psi(js), ygeo(js), yellip(js),
     1            yinden(js), ytrian(js), yshift(js)
            END IF

         END DO
      END DO

   95 FORMAT(/,71('-'),/,' Toroidal Plane: Phi = 180/Nfp',/,71('-'),/)
  120 FORMAT(1x,i5,6f14.5,1p,3e14.2)

      WRITE (nthreed, 130)
  130 FORMAT(//,' Magnetic Fields and Pressure',/,1x,71('-'))
      fac = p5/mu0
      WRITE (nthreed, 140) sump/mu0, pavg/mu0, fac*sumbpol, 
     1   fac*sumbpol/volume_p, fac*sumbtor, fac*sumbtor/volume_p, 
     2   fac*sumbtot, fac*sumbtot/volume_p, c1p5*sump/mu0, 
     3   c1p5*pavg/mu0
  140 FORMAT(' Volume Integrals (Joules) and Volume ',
     1   'Averages (Pascals)',/,24x,'Integral',6x,'Average',/,
     2   ' pressure         = ',1p,2e14.6,/,' bpol**2 /(2 mu0) = ',
     3   2e14.6,/,' btor**2/(2 mu0)  = ',2e14.6,/,
     4   ' b**2/(2 mu0)     = ',2e14.6,/,' EKIN (3/2p)      = ',
     5   2e14.6,/)

      WRITE (nthreed, 800)
  800 FORMAT(/,' MAGNETIC AXIS COEFFICIENTS'/,
     1   '    n     rmag       zmag        rmag        zmag',/,
     2   '        (cos nv)   (sin nv)    (sin nv)    (cos nv)',/)
      loff = LBOUND(rmags,1)
      DO n = 0, ntor
         n1 = n + loff
         t1 = mscale(0)*nscale(n)
         tz = t1
         IF (.not.lthreed) tz = 0
         IF (lasym) THEN
            WRITE (nthreed, 820) n, t1*rmags(n1), (-tz*zmags(n1)),
     1                             -tz*rmaga(n1),   t1*zmaga(n1)
         ELSE
            WRITE (nthreed, 820) n, t1*rmags(n1), (-tz*zmags(n1))
         END IF
      END DO
  820 FORMAT(i5,1p,4e12.4)

      IF (lrecon) THEN

      betstr = two*SQRT(sump2/volume_p)/(sumbtot/volume_p)

      WRITE (nthreed, 150) betatot, betapol, betator
  150 FORMAT(' From volume averages over plasma, betas are',/,
     1   ' beta total    = ',f14.6,/,' beta poloidal = ',f14.6,/,
     2   ' beta toroidal = ',f14.6,/)

      WRITE (nthreed, 160) 1.e-6_dp*ctor/mu0, rbtor, betaxis, betstr
  160 FORMAT(' Toroidal Current     = ',f14.6,'  [MA]',/
     1   ' R * Btor-vac         = ',f14.6,' [Wb/M]',/,
     2   ' Peak Beta            = ',f14.6,/,' Beta-star            = ',
     3   f14.6,/)
!
!
!     Shafranov surface integrals s1,s2
!     Plasma Physics vol 13, pp 757-762 (1971)
!     Also, s3 = .5*S3, defined in Lao, Nucl. Fusion 25, p.1421 (1985)
!     Note: IF ctor = 0, USE Int(Bsupu*Bsubu dV) for ctor*ctor/R
!
      IF (lfreeb) THEN
        factor = zero                             !Compute current-like norm
        DO l = ns, nrzt, ns
           js = MOD(l - 1,ns) + 1
           ncol = (l - 1)/ns + 1
           factor = factor + twopi*wint(l)*ABS(bsubu(js,ncol,0))
        END DO
        factor = one/factor**2
        facnorm = factor*twopi*twopi

        ALLOCATE (redg1u(nznt), rbps1u(nznt))
        redg1u(:nznt) = r1(ns:nznt*ns:ns,0) + r1(ns:nznt*ns:ns,1)
        rbps1u(:nznt) = two*facnorm*redg1u(:nznt)*(bpolvac(:nznt) +
     1    presf(ns))*wint(ns:nznt*ns:ns)
        sigr0 = DOT_PRODUCT(rbps1u(:nznt),zu0(ns:nznt*ns:ns))
        sigr1 = DOT_PRODUCT(rbps1u(:nznt)*zu0(ns:nznt*ns:ns),
     1                    redg1u(:nznt))
        sigz1 = -DOT_PRODUCT(rbps1u(:nznt)*ru0(ns:nznt*ns:ns),
     1           z1(ns:nznt*ns:ns,0) + z1(ns:nznt*ns:ns,1))
        DEALLOCATE (redg1u, rbps1u)

        er = sigr1 + sigz1
        rlao = volume_p/(twopi*cross_area_p)               !LAO, NUCL.FUS.25(1985)1421
        flao = rshaf/rlao
        fgeo = rshaf/rcen
        factor = two*factor/rshaf

        smaleli = factor*sumbpol
        betai = two*factor*sump
        musubi = vnorm*factor*musubi
        dmusubi_meas = two*twopi*factor*delphid*rbtor
        lambda = p5*smaleli + betai
        s11 = (er - rshaf*sigr0)/rshaf         !Shafranov def. based on RT
        s12 = (er - rcen*sigr0)/rcen               !R = Rgeometric
        s13 = (er - rlao*sigr0)/rlao               !R = RLao
        s2 = sigr0
        s3 = sigz1/rshaf
        delta1 = zero
        delta2 = one - fgeo
        delta3 = one - flao
        WRITE (nthreed, 170) rshaf, rcen, rlao, dmusubi_meas, delta1,
     1   delta2, delta3, s11, s12, s13, s2, s2, s2, s3, s3*fgeo, s3*flao
     2   , smaleli, smaleli*fgeo, smaleli*flao, musubi, musubi*fgeo,
     3   musubi*flao, betai, betai*fgeo, betai*flao, musubi + s11,
     4   musubi*fgeo + s12 + s2*(one - fgeo), musubi*flao + s13 + s2*(
     5   one - flao), lambda, fgeo*dlam, flao*dlam, 0.5*s11 + s2, 0.5*(
     6   s12 + s2*(one + fgeo)), 0.5*(s13 + s2*(one + flao)), (3*betai
     7    + smaleli - musubi)*0.5/(s11 + s2) - one, fgeo*(3*betai +
     8   smaleli - musubi)*0.5/(s12 + s2) - one, flao*(3*betai + smaleli
     9    - musubi)*0.5/(s13 + s2) - one, (betai + smaleli + musubi)*0.5
     .   /s2 - one, fgeo*(betai + smaleli + musubi)*0.5/s2 - one,
     .   flao*(betai + smaleli + musubi)*0.5/s2 - one
  170 FORMAT(' Integrals of Shafranov',/,1x,22('-'),/,
     1   ' RT (Flux-weighted)   = ',f14.6,' [M]',/,
     2   ' RG (Geometric)       = ',f14.6,' [M]',/,
     3   ' RL (Vol/2*pi*Area)   = ',f14.6,' [M]',/,
     4   ' Mui (diamagnetism)   = ',f14.6,2/,32x,'R = RT',12x,'R = RG',
     5   12x,'R = RL',/,20x,3(10x,8('-')),/,' delta = 1 - RT/R     = ',3
     6   (f14.6,4x),/,' s1                   = ',3(f14.6,4x),/,
     7   ' s2                   = ',3(f14.6,4x),/,
     8   ' s3                   = ',3(f14.6,4x),/,
     9   ' Li                   = ',3(f14.6,4x),/,
     .   ' Mui                  = ',3(f14.6,4x),/,
     1   ' Betai (Calculated)   = ',3(f14.6,4x),/,
     2   ' Betai (Mui + s1)     = ',3(f14.6,4x),/,
     3   ' Lambda (Calculated)  = ',3(f14.6,4x),/,
     4   ' Lambda (s1/2 + s2)   = ',3(f14.6,4x),/,
     5   ' 1st Shafr''v relation = ',3(f14.6,4x),/,
     6   ' (3*Betai + Li - Mui)/[2*(s1+s2)] - 1',/,
     7   ' Radial force balance = ',3(f14.6,4x),/,
     8   ' (Betai + Li + Mui)/(2*s2) - 1',/)

      END IF
!
!     Safety Factors (q)
!
      qzero = one/ABS(iotaf(1))
      qedge = one/ABS(iotaf(ns))
      WRITE (nthreed, 180) qzero, qedge, qzero/qedge
  180 FORMAT(' Safety Factors',/,1x,14('-'),/,' q (on axis) = ',f12.4,/,
     1   ' q (at edge) = ',f12.4,/,' q(0)/qedge  = ',f12.4,/)

!
!     PRINT OUT IOTA, PRESSURE SPLINE COEFFICIENTS
!     (PRESSURE IN MKS UNITS, NWT/M**2)
!
      WRITE (nthreed, 190)
  190 FORMAT(/,' SUMMARY OF IOTA AND PRESSURE SPLINES'/,
     1   '  K   Spline Node       IOTA(K)      IOTA"(K)'/,
     2   '        SQRT(s)',/,3('-'),3(4x,10('-')))
      DO i = 1, isnodes
         WRITE (nthreed, 200) i, sknots(i), ystark(i), y2stark(i)
      END DO
  200 FORMAT(i3,1p,3e14.3)
      WRITE (nthreed, 210)
  210 FORMAT(/,'  K   Spline Node       PRES(K)      PRES"(K)'/,
     1   '        SQRT(s)',/,3('-'),3(4x,10('-')))
      factor = pthommax
      DO i = 1, ipnodes
         WRITE (nthreed, 200) i, pknots(i), factor*ythom(i), factor*
     1      y2thom(i)
      END DO

!
!     PRINT-OUT MAGNETIC PITCH DATA
!
      nchimse = imse
      total_mse_chi = zero

      IF (imse .ne. 0) THEN

         WRITE (nthreed, 220)
  220    FORMAT(//,4x,'N     R(data)      Spline',3x,
     1      'atan[BZ/BT] ATAN[BZ/BT] Chi-Sq-Err',4x,
     2      'Sigma       BZ(in)     BZ(calc)       BT        q(in)',/,
     3      12x,'[m]',5x,'SQRT(s)-node',2x,'(in-deg)',2x,'(calc-deg)',
     4      20x,3(9x,'[T]')/,2x,3('-'),1x,2(3x,9('-')),8(2x,10('-'))/)

         msewgt = zero
         denwgt = zero
         DO n = 1, imse + 1
            nsort = isortr(n)
            pitchc = ATAN(starkcal(n))/dcon
            pitchm = ATAN(datastark(nsort))/dcon
            wght0 = ATAN(sigma_stark(nsort))/dcon
            js = indexs1(n)
            lt = indexu1(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            rcalc = (one - delso1(n))*rjs + delso1(n)*rjs1
            IF (delse1(n) .lt. zero) rcalc = zero
            IF (rcalc .ne. 0.) btor = fpsical(n)/rcalc
            bzin = TAN(dcon*pitchm)*btor
            bzcalc = TAN(dcon*pitchc)*btor
            chisq = zero
            IF (ABS(rcalc - router) .lt. epstan) THEN
               WRITE (nthreed, 230) n, rcalc, rsort0(n), pitchm,
     1            pitchc, wght0, one/(qcalc(n) + epstan)
            ELSE
               IF (ABS(rsort(n) - rcalc) .le. epstan) THEN
                  chisq = ((pitchc - pitchm)/wght0)**2
                  msewgt = msewgt + chisq
                  denwgt = denwgt + (pitchm/wght0)**2
               ENDIF
               WRITE (nthreed, 240) n, rcalc, rsort0(n), pitchm,
     1            pitchc, chisq, wght0, bzin, bzcalc, btor,
     2            one/(qmeas(n) + epstan)
            ENDIF
         END DO

         chiwgt = msewgt/(imse)
         msewgt = SQRT(msewgt/denwgt)
c                                             !total chi-squared for mse
         total_mse_chi = (imse)*chiwgt
         WRITE (nthreed, 250) chiwgt, msewgt

  230    FORMAT(i5,' ',1p,2e12.3,2e12.3,12x,e12.3,6x,
     1      '- Outer (phantom) Edge -',6x,e12.3)
  240    FORMAT(i5,' ',1p,2e12.3,8e12.3)
  250    FORMAT(/' MEAN CHI-SQ ERROR IN STARK DATA MATCH : ',1p,e10.3,/,
     1      ' RMS ERROR IN STARK DATA MATCH : ',e10.3,2/)

      ENDIF

!
!     PRINT-OUT PRESSURE DATA
!
      nchipres = itse
      total_pres_chi = zero

      IF (lpprof) THEN
         tswgt = zero
         denwgt = zero
         WRITE (nthreed, 300) presfac*pfac
  300    FORMAT(4x,'N     R(in)      R(calc)',4x,f6.2,
     1      ' X Pres(in)      Pres(calc)','  Chi-Sq Err       Sigma'/,2x
     2      ,3('-'),2(2x,10('-')),2(6x,12('-')),2(3x,9('-')),/)

         DO n = 1, itse
            js = indexs2(n)
            lt = indexu2(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            IF (delso2(n) .eq. (-one)) delso2(n) = one
            rcalc = (one - delso2(n))*rjs + delso2(n)*rjs1
            wght0 = sigma_thom(n)
            chisq = ((datathom(n)*pfac-pcalc(n)/mu0)/wght0)**2
            tswgt = tswgt + chisq
            denwgt = denwgt + (datathom(n)/wght0)**2
            WRITE (nthreed, 310) n, rthom(n), rcalc, presfac*pfac*
     1         datathom(n), pcalc(n)/mu0, chisq, wght0
         END DO

         chiwgt = tswgt/(itse)
         tswgt = SQRT(tswgt/denwgt)
         total_pres_chi = (itse)*chiwgt !total
         WRITE (nthreed, 320) chiwgt, tswgt
  310    FORMAT(i5,1p,2e12.3,2e18.3,2e12.3)
  320    FORMAT(/' MEAN CHI-SQ ERROR IN PRESSURE DATA MATCH: ',
     1      1p,e10.3,/,' RMS ERROR IN PRESSURE DATA MATCH: ',e10.3/)
      ENDIF

!
!     SUMMARIZE MAGNETICS DATA AND MATCH
!
      CALL magnetics_data

!
!     COMPUTE REAL TOROIDAL CURRENT ALONG MIDPLANE (MULTIPLY IZETA BY R/SQRT(G))
!     IN PHI = 0 PLANE
!
      CALL getcurmid (curmid, izeta, gsqrt, r12)

      DO nout = nthreed, nmac, (nmac - nthreed)
         IF (nout .eq. nmac) THEN
            IF (.not.lmac) CYCLE
            WRITE (nout, *)
     1      'FOLLOWING DATA EQUALLY SPACED IN TOROIDAL FLUX'
         END IF
         WRITE (nout, 700)
         IF (nout .eq. nthreed) WRITE (nout, 710)
         iotas(1) = two*iotas(2) - iotas(3)
         iotas(ns+1) = two*iotas(ns) - iotas(ns1)
         vp(1) = two*vp(2) - vp(3)
         vp(ns+1) = two*vp(ns) - vp(ns1)
         pres(1) = two*pres(2) - pres(3)
         pres(ns+1) = two*pres(ns) - pres(ns1)
         DO icount = 1, 2*ns - 1
            js = MOD(imid(icount) - 1,ns) + 1
            ageo(icount) = ygeo(js)
            phimid(icount) = toroidal_flux*(js - 1)/(ns1)
            psimid(icount) = psi(js)
            volpsi(icount) = vnorm*SUM(vp(2:js))
            IF (js .eq. 1) volpsi(icount) = zero
            dlogidsi = (iotas(js+1)-iotas(js))*ohs/iotaf(js)
            vprime = p5*vnorm*ohs*(vp(js)+vp(js+1))
            pprime = (pres(js+1)-pres(js))*ohs
            rgeo = rmid(icount) + ygeo(js)
            IF (icount .gt. ns) rgeo = rgeo - two*ygeo(js)
            alfa(icount) = -two*pprime*vprime*SQRT(two*volpsi(icount)/
     1         rgeo/twopi)/(iotaf(js)*phipf(js))**2
            shear(icount) = -two*volpsi(icount)*dlogidsi/vprime
            WRITE (nout, 720) rmid(icount), ygeo(js), psi(js),
     1         volpsi(icount), qmid(icount), shear(icount),
     2         presmid(icount),
     3         alfa(icount), curmid(icount), datamse(icount)
         END DO
      END DO
  700 FORMAT(/3x,'RMID[M]',6x,'AMID[M]',6x,'PSI[Wb]',4x,'VOL[M**3]',9x,
     1   'q(r)',2x,'SHEAR(Spsi)',6x,'P(PASC)',6x,'ALF(P'')',2x,
     2   'JTOR(A/M**2)',1x,'ATAN(Bz/Bt)[deg]')
  710 FORMAT(10('-'),9(3x,10('-')))
  720 FORMAT(1p,e10.3,9(3x,e10.3))

!     Determine q-MIN on the s grid by MIN-splining on the SQRT-s grid
      nptsx = 2*ns - 1

      wmidx(:nptsx) = one
      tenmidx(:nptsx) = 0.1_dp
      rmidx(:nptsx) = rmid(:nptsx)
      qmidx(:nptsx) = qmid(:nptsx)

      tol = .005_dp
      sqmin = fmax(zero,one,splints,tol)
      iotamax = splints(sqmin)
      qmin3 = -99999.0_dp
      IF (iotamax .ne. zero) qmin3 = one/iotamax
      sqmin = sqmin**2

!     Determine q-MIN on a fine r-midplane grid
      tol = .005_dp
!     outboard side
      rqmin1 = fmin(rmid(ns),rmid(nptsx),splintx,tol)
      qmin1 = splintx(rqmin1)
      rqmin2 = fmin(rmid(1),rmid(ns),splintx,tol)!outboard side only
      qmin2 = splintx(rqmin2)
      WRITE (nthreed, 730) qmin1, rqmin1, qmin2, rqmin2, qmin3, sqmin
  730 FORMAT(//' MINIMUM Q :     OUTBOARD SIDE: QMIN= ',f6.3,'  AT R= ',
     1   f6.3,/,'                  INBOARD SIDE: QMIN= ',f6.3,'  AT R= '
     2   ,f6.3,/,'                    IN S SPACE: QMIN= ',f6.3,
     3   '  AT S= ',f6.3,/)

      nchiiota0 = 0
      total_chi_square = total_b_chi + total_saddle_chi + total_pres_chi
     1    + total_mse_chi + total_chi_cur + total_chi_dia + chisqerr(
     2   islope0)                        !Total CHISQ of equilibrium fit
      nchitot = nbfldn + nchiiota0 + nchisaddle + nchipres + nchimse +
     1   nchicur + nchidia
      total_chi_square_n = total_chi_square/MAX(1,nchitot)
      WRITE (nthreed, 900)
      WRITE (nthreed, 901) (bloopnames(n),nbfld(n),b_chi(n),b_chi(n)/
     1   MAX(1,nbfld(n)),n=1,nbsets)
      WRITE (nthreed, 902) nchisaddle, total_saddle_chi,
     1   total_saddle_chi/MAX(1,nchisaddle), nchipres, total_pres_chi,
     2   total_pres_chi/MAX(1,nchipres), nchimse, total_mse_chi,
     3   total_mse_chi/MAX(1,nchimse), nchicur, total_chi_cur,
     4   total_chi_cur/MAX(1,nchicur), nchidia, total_chi_dia,
     5   total_chi_dia, nchitot, total_chi_square, total_chi_square_n

  900 FORMAT(/,' CHI-SQUARED SUMMARY:',t25,'Type',t50,'Nchi',t65,'ChiSq'
     1   ,t84,'ChiSq/N'/,t25,'----',t50,'----',t65,'-----',t84,'-------'
     2   )
  901 FORMAT(t25,'B-Loops-',a,t50,i3,t60,1p,e12.4,t80,e12.4)
  902 FORMAT(t25,'Saddle',t50,i3,t60,1p,e12.4,t80,e12.4,/,t25,
     1   'Pressure',t50,i3,t60,e12.4,t80,e12.4,/,t25,'MSE',t50,i3,
     2    t60,e12.4,t80,e12.4,/,t25,'Ip',t50,i3,t60,e12.4,t80,
     3    e12.4,/,t25,'Diamagnetic Flux',t50,i3,t60,e12.4,t80,e12.4,/
     4   ,t25,'TOTAL',t50,i3,t60,e12.4,t80,e12.4)

      END IF          !!IF(LRECON)


      END SUBROUTINE eqfor
