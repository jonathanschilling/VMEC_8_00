      SUBROUTINE bcovar (lu, lv, tau, lmnsc, lscreen)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp
      USE realspace
      USE vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1             rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2             bsubu_e => clmn_e, bsubv_e => blmn_e, 
     3             bsubu_o => clmn_o, bsubv_o => blmn_o,
     4             bsq => bzmn_o, phipog => brmn_o
      USE vsvd, ONLY: phifac, phifsave, imovephi
      USE xstuff, ONLY: xc
      USE precon2d, ONLY: ljog_test, loglam_blks, lprec2d
      USE fbal
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt,0:1), INTENT(inout) :: lu, lv
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), 
     1     INTENT(inout) :: lmnsc
      REAL(rprec) :: tau(nrzt)
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(rprec), PARAMETER :: c1p5 = (one + p5)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, is, js, ndim
      REAL(rprec) :: r2, r3
      REAL(rprec) :: arnorm, aznorm, volume, tcon_mul
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: luu, luv, lvv
      REAL(rprec) :: curtor_temp, curpol_temp
      REAL(rprec), DIMENSION(:), POINTER :: bsupu, bsubuh, 
     1                                      bsupv, bsubvh, r12sq
      LOGICAL :: lflam
C-----------------------------------------------
      ndim = 1+nrzt
      lflam = (iequi.ne.2 .or. ljog_test .or. lasym)   

!
!     POINTER ALIAS ASSIGNMENTS
!   
      bsupu => bsubu_o;  bsubuh => bsubu_o
      bsupv => bsubv_o;  bsubvh => bsubv_o
      r12sq => bsq

      ALLOCATE (luu(ndim), luv(ndim), lvv(ndim), stat=is)

      IF (is .ne. 0) STOP 'allocation error in bcovar'

!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      guv = 0;  gvv = 0

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
!

CDIR$ IVDEP
      DO l = 1, nrzt
         r2 = sqrts(l)*sqrts(l)
         guu(l)   =  ru(l,0)*ru(l,0) + r2*ru(l,1)*ru(l,1)
     1            +  zu(l,0)*zu(l,0) + r2*zu(l,1)*zu(l,1)
         luu(l)   = (ru(l,0)*ru(l,1) + zu(l,0)*zu(l,1))*2
         r12sq(l) = r1(l,0)*r1(l,0) + r2*r1(l,1)*r1(l,1)              !Comment: r12sq = r12**2
         phipog(l)= 2*r1(l,0)*r1(l,1)                                 !Comment: r12sq = r12**2
      END DO

      IF (lthreed) THEN
CDIR$ IVDEP
         DO l = 1, nrzt
            r2 = sqrts(l)*sqrts(l)
            guv(l)   =  ru(l,0)*rv(l,0) + r2*ru(l,1)*rv(l,1)
     1               +  zu(l,0)*zv(l,0) + r2*zu(l,1)*zv(l,1)
            luv(l)   =  ru(l,0)*rv(l,1) + ru(l,1)*rv(l,0)
     1               +  zu(l,0)*zv(l,1) + zu(l,1)*zv(l,0)
            gvv(l)   =  rv(l,0)*rv(l,0) + r2*rv(l,1)*rv(l,1)
     1               +  zv(l,0)*zv(l,0) + r2*zv(l,1)*zv(l,1)
            lvv(l)   = (rv(l,0)*rv(l,1) + zv(l,0)*zv(l,1))*2
         END DO
      END IF
      DO l = nrzt, 2, -1
         guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
         r12sq(l) = p5*(r12sq(l) + r12sq(l-1) + shalf(l)*             !Comment: r12sq = r12**2
     1                (phipog(l) + phipog(l-1)))                      !Comment: r12sq = r12**2
      END DO
      IF (lthreed) THEN
         DO l = nrzt, 2, -1
            guv(l) = p5*(guv(l) + guv(l-1) +
     1         shalf(l)*(luv(l) + luv(l-1)))
            gvv(l) = p5*(gvv(l) + gvv(l-1) +
     1         shalf(l)*(lvv(l) + lvv(l-1)))
         END DO
      END IF

      tau(1:nrzt) = gsqrt(1:nrzt)
      gsqrt(1:nrzt) = r12(1:nrzt)*tau(1:nrzt)      

      gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)

      DO l = 1, nrzt
         IF (gsqrt(l) .eq. zero) THEN
            phipog(l) = 0
         ELSE
            phipog(l) = phip(l)/gsqrt(l)
         END IF
      END DO

      phipog(1)    = 0
      phipog(ndim) = 0

      DO js = 2, ns
         vp(js) = signgs*SUM(gsqrt(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO
      IF (iter2 .eq. 1) voli = twopi*twopi*hs*SUM(vp(2:ns))

!
!     RECONSTRUCT IOTA, PRESSURE PROFILE FROM DATA
!
      IF (lrecon) CALL newprofil (phipog)
      IF (.not.lrecon .and. imovephi.gt.0) CALL newphi (phipog)

!
!     STORE CURRENT(0), MAGNETIC PITCH, PRESSURE DATA
!     STORE HALF-MESH VALUES OF LU, LV FOR ACCURATE CURRENT CALCULATION
!
      IF (iequi .eq. 1) CALL storesvd (r1(1,0), r1(1,1), lu(1,0),
     1   lu(1,1), lv(1,0), phipog, zu0)

!
!     NOTE: LU = 1+LAMU, LV = -LAMV COMING INTO THIS ROUTINE
!     WILL ADD IOTAF AFTER CALL TO GETIOTA...
!          BSUPU = (PHIP/GSQRT)*(iota - LAMV),
!          BSUPV = (PHIP/GSQRT)*(1    + LAMU)
!

!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!
CDIR$ IVDEP
      DO l = 2, nrzt
         bsupu(l) = p5*phipog(l)*(lv(l,0) + lv(l-1,0) + shalf(l)*
     1                           (lv(l,1) + lv(l-1,1)))
         bsupv(l) = p5*phipog(l)*(lu(l,0) + lu(l-1,0) + shalf(l)*
     1                           (lu(l,1) + lu(l-1,1)))
      END DO

      bsupv(1) = 0
      bsupu(1) = 0

!
!     UPDATE IOTA EITHER BY SOLVING <Bsubu> = icurv EQUATION
!     OR BY PUSHING IOTA IN TIME, USING Force-iota = <Bsubu> - icurv 
      IF (.not.liota) THEN
!
!     COMPUTE UPDATED IOTA PROFILE
!
         CALL getiota(phipog, bsupu, bsupv)

!
!        STORE iotaf IN lmnsc(m=0,n=0) TO BE EVOLVED WHEN 2d PRECONDITIONER IS USED
!        EXPERIENCE SHOWS liota = .false. WORKS BEST; LEFT IT IN JUST FOR FUTURE
!        POSSIBLE USE. MODIFICATIONS TO ANALYTIC HESSIAN FOR LFORBAL = TRUE ARE NOT IMPLEMENTED
!
!         liota = (ncurr.eq.1) .and. ((iter2.gt.50 .and. fsq.lt.1.E-6_dp)
!     1                        .or. lprec2d)   !TURN ON FOR TESTING....
         IF (liota) lmnsc(1:ns,0,0) = iotaf(1:ns)      
         IF (liota) PRINT *,' Turning on lambda force'
     
      ELSE IF (ncurr .eq. 1) THEN
         iotaf(2:ns) = lmnsc(2:ns,0,0)
         iotaf(1) = 2*iotaf(2) - iotaf(3)     !Note: this breaks tri-diagonal coupling in Hessian 
         iotas(2:ns) = (iotaf(2:ns) + iotaf(1:ns-1))/2
      END IF

!     ADD PRESENT VALUE OF IOTA HERE.
      DO js = 1, ns
         bsupu(js:nrzt:ns) = bsupu(js:nrzt:ns) 
     1                     + iotas(js)*phipog(js:nrzt:ns)
      END DO

!
!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY
!
      IF (lflam) THEN
         lvv = phipog(:ndim)*gvv
         bsubv_e(1:nrzt) = p5*(lvv(1:nrzt)+lvv(2:nrzt+1))*lu(1:nrzt,0)
      END IF

!     STORE CONTRAVARIANT B COMPONENTS IN LU, LV
      lu(:nrzt,0) = bsupv(:nrzt)
      lv(:nrzt,0) = bsupu(:nrzt)

      bsubvh(:nrzt) = guv(:nrzt)*lv(:nrzt,0)

      IF (lflam) THEN
         lvv = lvv*shalf
         bsubv_e(1:nrzt) = bsubv_e(1:nrzt) + p5*
     1                    ((lvv(1:nrzt) + lvv(2:nrzt+1))*lu(1:nrzt,1)
     2                   + (bsubvh(1:nrzt) + bsubvh(2:nrzt+1)))
      END IF

!
!     FINALLY, COMPUTE LAMBDA FORCES (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
      bsubuh(:nrzt) = guu(:nrzt)*lv(:nrzt,0) + guv(:nrzt)*lu(:nrzt,0)
      bsubvh(:nrzt) = bsubvh(:nrzt) + gvv(:nrzt)*lu(:nrzt,0)

      bsubuh(ndim) = 0
      bsubvh(ndim) = 0

!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
      bsq(:nrzt) = p5*(lv(:nrzt,0)*bsubuh(:nrzt) +
     1                 lu(:nrzt,0)*bsubvh(:nrzt))
      wb = hs*ABS(SUM(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*SUM(vp(2:ns)*pres(2:ns))

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
      CALL calc_fbal(bsubuh, bsubvh)
    
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      ctor = signgs*twopi*SUM((c1p5*bsubuh(ns:nrzt:ns) -
     1              p5*bsubuh(ns-1:nrzt:ns))*wint(ns:nrzt:ns))

!     TESTING....
      IF (((MOD(iter2,nstep).eq.0 .and. iequi.eq.0) .or. iequi.eq.1)
     1      .and. lscreen) THEN
        WRITE (6,'(2x,a,1p,e12.3)')
     1         'Force balance: <(JXB - grad_p)> = ', 
     2         SUM(ABS(equif(2:ns-1)/fpsi(2:ns-1))) * hs**2
      END IF
!     END TESTING...

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!

      IF (lflam) THEN
         DO l = 1, nrzt
            r2 = (1 - sqrts(l)*sqrts(l))**2 * pdamp 
            r3 = p5*(1 - r2)
            bsubu_e(l) = p5*(bsubuh(l) + bsubuh(l+1))    
            bsubv_e(l) = bsubv_e(l)*r2
     1                 + r3*(bsubvh(l) + bsubvh(l+1))
         END DO
      END IF

!
!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
!
      DO js = 2,ns
         bsq(js:nrzt:ns) = bsq(js:nrzt:ns) + pres(js)
      END DO

!TESTING....
      IF (MOD(iter2,nstep).eq.0 .and. liota) THEN
         r3 = 0
         DO js = 2, ns
            r2 =  SUM(wint(js:nrzt:ns)*bsubu_e(js:nrzt:ns)) 
            r3 = r3 + (r2 - p5*(icurv(js) + icurv(js-1)))**2
         END DO
         PRINT *,' <F-iota>**2 = ', r3
      END IF
!END TESTING
!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS
!
!      IF (iequi.eq.2 .and. loglam_blks) 
!     1   CALL lam_blks (phipog, guu, guv, gvv, tau)

      IF (MOD(iter2-iter1,ns4).eq.0 .and. iequi.eq.0 
     1        .and. .not.lprec2d) THEN
         phifsave = phifac
         phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt)
         CALL lamcal(phipog, guu, guv, gvv)
         CALL precondn(lu,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,crd,rzu_fac,cos01)
         CALL precondn(lu,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                 r1(1,1),azm,azd,bzm,bzd,crd,rru_fac,sin01)

         rzu_fac(2:ns-1) = sqrts(2:ns-1)*rzu_fac(2:ns-1)
         rru_fac(2:ns-1) = sqrts(2:ns-1)*rru_fac(2:ns-1)
         frcc_fac(2:ns-1) = one/rzu_fac(2:ns-1);  rzu_fac = rzu_fac/2
         fzsc_fac(2:ns-1) =-one/rru_fac(2:ns-1);  rru_fac = rru_fac/2

         guu(:ndim) = guu(:ndim)*r12(:ndim)**2
         volume = hs*SUM(vp(2:ns))
         r2 = MAX(wb,wp)/volume
         fnorm = one/(SUM(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))
         fnorm1 = one/SUM(xc(1+ns:2*irzloff)**2)

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 value from old-style file
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((2*r0scale**2)*(4*r0scale**2))       !!Scaling of ard, azd (2*r0scale**2); 
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)

         DO js = 2, ns-1
           arnorm = SUM(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
           aznorm = SUM(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
           IF (arnorm.eq.zero .or. aznorm.eq.zero)
     1        STOP 'arnorm or aznorm=0 in bcovar'

           tcon(js) = MIN(ABS(ard(js,1)/arnorm),
     1                    ABS(azd(js,1)/aznorm)) * tcon_mul*(32*hs)**2
         END DO
         tcon(ns) = p5*tcon(ns-1)
      ENDIF

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      IF (iequi .eq. 1) THEN

!DBG     luu = bsubuh
!DBG     lvv = bsubvh

         DO js = ns-1,2,-1
            DO l = js, nrzt, ns
!blend              bsubuh(l) = 2*bsubu_e(l) - bsubuh(l+1)
               bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
            END DO
         END DO

!     ADJUST <bsubvh> AFTER MESH-BLENDING
         DO js = 2, ns
!blend             curtor_temp = jtor(js) 
!blend    1                  - SUM(bsubuh(js:nrzt:ns)*wint(js:nrzt:ns))
            curpol_temp = fpsi(js) 
     1                  - SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))
            DO l = js, nrzt, ns
!blend             bsubuh(l) = bsubuh(l) + curtor_temp
               bsubvh(l) = bsubvh(l) + curpol_temp
            END DO
         END DO

!DBG         DO js = 2, ns
!DBG            curtor_temp = SUM(bsubuh(js:nrzt:ns)*wint(js:nrzt:ns))
!DBG            curpol_temp = SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))
!DBG            WRITE (33,'(a,i4,2(a,1pe12.3),/,2(a,1pe12.3))')
!DBG     1     ' js = ',js,' CURTOR = ', curtor_temp,' JCURV = ',icurv(js),
!DBG     1              ' CURPOL = ', curpol_temp, ' FPSI = ', fpsi(js)
!DBG         END DO

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)

!DBG     DO js = 2, ns
!DBG        PRINT *, 'JS = ', js
!DBG        PRINT *,
!DBG 1      '    BSUBU(old)     BSUBU(new)    BSUBV(old)     BSUBV(new)'
!DBG        DO l = js, nrzt, ns
!DBG           WRITE(*,1223) luu(l), bsubu_e(l), lvv(l), bsubv_e(l)
!DBG        END DO
!DBG     END DO

 1223    FORMAT(1p,4e14.5)

         GOTO 1000

      END IF

!
!     MAKE HESSIAN APPROXIMATELY SYMMETRIC AND POSITIVE DEFINITE
!     (NOTE: DIVIDE BZ1 BY PHIP(2)**2 IN RESIDUE, AND MULT FACLAM BY PHIP(2))
!
      IF (lflam) THEN
         bsubu_e(1:nrzt) =-phip(2)*bsubu_e(1:nrzt)
         bsubv_e(1:nrzt) =-phip(2)*bsubv_e(1:nrzt)

         bsubu_o(:nrzt) = sqrts(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = sqrts(:nrzt)*bsubv_e(:nrzt)
      END IF

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
      DO l=2,nrzt
         guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
         guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
         gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)
         lv(l,0)  = bsq(l)*tau(l)
         lu(l,0)  = bsq(l)*r12(l)
      ENDDO

 1000 CONTINUE

      DEALLOCATE (luu, luv, lvv, stat=l)

      END SUBROUTINE bcovar
