      SUBROUTINE funct3d (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac
      USE vmec_params, ONLY: bad_jacobian_flag
      USE realspace
      USE vforces
      USE vsvd, ONLY: router, rinner, gphifac, grmse
      USE xstuff
      USE timer_sub
      USE precon2d
      USE vmec_utils, ONLY: cyl2flx
      USE vparams, ONLY: twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l0pi, l, lk, ivacskip
      INTEGER :: nvskip0 = 0
      REAL(rprec), DIMENSION(mnmax) ::
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(nznt) :: rax, zax
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: extra1
     1   , extra2, extra3, extra4
      REAL(rprec), DIMENSION(:), POINTER :: lu, lv, gcon
      REAL(rprec) :: presf_ns, delr_mse, delt0

      INTEGER :: i, j, k, nsum, nfmin, nfe, iflag
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin, fmin_max, s_in, u_in
      REAL(rprec) :: ton, toff, ttotal
C-----------------------------------------------
!
!     POINTER ALIASES
!
      lu => czmn;  lv => crmn;  gcon => z1(:,0)

      CALL second0 (tfunon)

      ALLOCATE (extra1(nrzt,0:1), stat=l)
      IF (lasym) ALLOCATE (
     1   extra2(nrzt,0:1), extra3(nrzt,0:1), extra4(nrzt,0:1),
     2   stat=l)
      IF (l .ne. 0) STOP 'Allocation error in funct3d'
!
!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!
      gc(:neqs2) = xc(:neqs2)*scalxc(:neqs2)

!
!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!
      IF (lrecon) THEN
         delr_mse = xc(neqs2)
         gc(1:ns) = gc(1:ns) + delr_mse
      ENDIF

!
!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!
      CALL second0 (tffton)
      CALL totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

!
!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
      IF (lasym) THEN

         CALL totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4, 
     1                 blmn, clmn, extra1, extra2)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, CONSTRAINTS ON FULL RANGE OF u (0 to 2*pi)

         CALL symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon, armn,
     1   brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2)

      ENDIF


      CALL second0 (tfftoff)
      IF (iequi .lt. 2) timer(tfft) = timer(tfft) + (tfftoff - tffton)

!
!     IN HESSIAN LOOP, ONLY COMPUTE PERTURBATIONS TO R, Z FOR ONE m,n,ntype
!
      IF (iequi .eq. 2) THEN
         xcdot(:neqs2) = xcdot(:neqs2)*scalxc(:neqs2)    
         CALL totzsps_hess (xcdot, r1, ru, rv, z1, zu, zv, lu, lv, 
     1                      rcon, zcon)
      END IF

      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
      router = r1(ns,0) + r1(ns,1)
      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = r1(1,0)
      z00 = z1(1,0)

!
!     COMPUTE CONSTRAINT RCON, ZCON
!
      DO l = 1,nrzt
         rcon(l,0) = rcon(l,0) + rcon(l,1)*sqrts(l)
         zcon(l,0) = zcon(l,0) + zcon(l,1)*sqrts(l)
         ru0(l) = ru(l,0) + ru(l,1)*sqrts(l)
         zu0(l) = zu(l,0) + zu(l,1)*sqrts(l)
      END DO

!
!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS SO THAT RESTART FOR FIXED BOUNDARY DOES NOT
!     HAVE A DISCONTINUITY DUE TO NEW RCON0....
!
      IF (iter2.eq.1 .and. .not.lfreeb) THEN
         DO l = 1, ns
            rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
            zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
         END DO
!        rcon0(:nrzt) = rcon(:nrzt,0)
!        zcon0(:nrzt) = zcon(:nrzt,0)
      ENDIF

!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
      CALL jacobian
      IF (irst.eq.2 .and. iequi.eq.0) GOTO 100

!
!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
!
      CALL second0 (tbcovon)
      CALL bcovar (lu, lv, extra1, xc(1+2*irzloff), lscreen)
      CALL second0 (tbcovoff)
      IF (iequi .lt. 2) timer(tbcov) = timer(tbcov) 
     1                               + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).

      IF (lfreeb .and. iter2.gt.1) THEN
         IF (fsqr + fsqz .le. 1.e-1_dp) ivac = ivac + 1
         IF (nvskip0 .eq. 0) nvskip0 = MAX(1, nvacskip)
         IF (ivac .ge. 0) THEN
            CALL second0 (tvacon)
            ivacskip = MOD(iter2 - iter1,nvacskip)
            IF (ivac .le. 2) ivacskip = 0
!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .eq. 0) THEN
               nvacskip = 1/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
               nvacskip = MAX (nvacskip, nvskip0)
            END IF

!           IF (lprec2d) ivacskip = MAX(1, ivacskip)

            DO lk = 1, nznt
               rax(lk) = r1(1 + ns*(lk - 1),0)
               zax(lk) = z1(1 + ns*(lk - 1),0)
            END DO
            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)
            CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn, rax, zax,
     1         ctor, rbtor, wint, ns, ivacskip, ivac, mnmax,
     2         ier_flag, lscreen)
            IF (ier_flag .ne. 0) RETURN
!
!           RESET FIRST TIME FOR SOFT START
!
            IF (ivac .eq. 1) THEN
               irst = 2;  delt0 = delt
               CALL restart_iter(delt0)
               irst = 1
            END IF

!
!           IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!           UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!           IF NON-VARIATIONAL FORCES ARE DESIRED
!
            presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
!RPRES      IF (iequi .ne. 1) presf_ns = 0

            lk = 0
            DO l = ns, nrzt, ns
               lk = lk + 1
               bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
               rbsq(lk) = (bsqvac(lk) + presf_ns)*ohs*(r1(l,0)+r1(l,1))
               dbsq(lk) = ABS(bsqvac(lk) + presf_ns - bsqsav(lk,3))
            END DO
            IF (ivac .eq. 1) THEN
               bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
               bsqsav(:nznt,2) = bsqvac(:nznt)
            ENDIF
            CALL second0 (tvacoff)
            IF (iequi .ne. 2) timer(tvac) = timer(tvac) 
     1                                    + (tvacoff - tvacon)
         ENDIF
      ENDIF

!
!     COMPUTE CONSTRAINT FORCE
!
      IF (iequi .ne. 1) THEN
         extra1(:nrzt,1) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     1                   + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
         CALL alias (gcon, extra1(:,0), extra1(:,1), gc, gc(1+mns),
     1               gc(1+2*mns))
      ELSE 
         IF (lrecon) xc(:ns) = xc(:ns) + delr_mse
         GOTO 100
      END IF

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      CALL second0 (tforon)
      CALL forces
!
!     SYMMETRIZE FORCES (in u-v space)
!
      IF (lasym) CALL symforce (armn, brmn, crmn, azmn, bzmn,
     1     czmn, blmn, clmn, rcon, zcon, r1, ru, rv, z1, zu, zv,
     2     extra3, extra4, extra1, extra2)

      CALL second0 (tforoff)
      IF (iequi .lt. 2) timer(tfor) = timer(tfor) + (tforoff - tforon)
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
      CALL second0 (tffton)
      CALL tomnsps (gc, armn, brmn, crmn, azmn, bzmn, czmn, 
     1              blmn, clmn, rcon, zcon)

      IF (lasym) CALL tomnspa (gc, r1, ru, rv, z1, zu, zv,
     1                         extra3, extra4, extra1, extra2)
      CALL second0 (tfftoff)
      IF (iequi .lt. 2) timer(tffi) = timer(tffi) + (tfftoff - tffton)

!================================================================
!
!     COMPUTE FORCE RESIDUALS. FOR IEQUI=2, THIS IS BEING CALLED
!     FROM COMPUTE_BLOCKS, AND MUST NOT BE PRECONDITIONED
!
!================================================================
      CALL second0 (treson)

      gc = gc * scalxc

      CALL residue (gc, gc(1+irzloff), gc(1+2*irzloff))

      CALL second0 (tresoff)
      IF (iequi .lt. 2) timer(tres) = timer(tres) + (tresoff - treson)

      gc(neqs1) = gphifac
      IF (iopt_raxis .gt. 0) gc(neqs2) = grmse

 100  CONTINUE

      DEALLOCATE (extra1)
      IF (lasym) DEALLOCATE (extra2, extra3, extra4)
      CALL second0 (tfunoff)
      IF (iequi .lt. 2) timer(tfun) = timer(tfun) + (tfunoff - tfunon)

      END SUBROUTINE funct3d
