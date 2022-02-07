      SUBROUTINE wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu,
     1   izeta, ier_flag, lwrite)
      USE vmec_main, p5 => cp5, two => c2p0
      USE vmec_params, ONLY: mscale, nscale, signgs, version_
      USE vmercier
      USE vmec_persistent
      USE vsvd
      USE vspline
      USE xstuff
      USE vmec_io
      USE realspace
      USE safe_open_mod
      USE mgrid_mod
      USE vacmod, ONLY: potvac,xmpot,xnpot,mnpd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu, izeta
      LOGICAL :: lwrite
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iasym, ireconstruct, i, j, js, lj, mn, lk, isgn,
     1   m, nmin0, n, k, iwout0, ntotal, n1, nwout, mn0
      REAL(rprec), DIMENSION(nznt) :: bmod
      REAL(rprec), DIMENSION(mnmax) :: rmnc, zmns, lmns,
     1   rmns, zmnc, lmnc, bmodmn, bmodmn1, lmnsjm1, lmncjm1
      REAL(rprec) :: dmult, gmn, bmn, currvmn, bsubumn, bsubvmn,
     1   bsubsmn, bsupumn, bsupvmn, tcosi, tsini, presfactor, 
     2   lmnsh, lmnch, gmns, bmns, bsubumns, bsubvmns, bsubsmnc,
     3   bsupumns, bsupvmns, potsin, potcos
      CHARACTER*120 :: wout_file,cdum
C-----------------------------------------------
      rmns = 0
      zmnc = 0
      lmnc = 0
!
!     THIS SUBROUTINE CREATES THE TEXT FILE WOUT WHICH
!     CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (FULL MESH), LMN (HALF MESH - CONVERTED FROM
!     INTERNAL FULL MESH REPRESENTATION)
!
!     IZETA (FULL), BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF)
!

      wout_file = 'wout_'//TRIM(input_extension)//'.txt'
      nwout = nwout0
      CALL safe_open(nwout, iwout0, wout_file, 'replace', 'formatted')
      IF (iwout0 .ne. 0) STOP 'Error writing WOUT file in VMEC WROUT'

      IF (lasym) THEN
         iasym = 1                                  ! asymmetric mode
      ELSE
         iasym = 0
      END IF
      IF (lrecon) THEN
         ireconstruct = 1
      ELSE
         itse = 0
         imse2 = 0
         ireconstruct = 0
      END IF

!
!     Insert version information into wout file. This will be parsed in
!     read_wout_file to RETURN the REAL value version_ to check the version number.
!
      WRITE (nwout, '(a15,a)') 'VMEC VERSION = ', version_
      WRITE (nwout, *) wb, wp, gamma, pfac,
     1  rmax_surf, rmin_surf, zmax_surf

      WRITE (nwout, *) nfp, ns, mpol, ntor, mnmax, itfsq, niter,
     1  iasym, ireconstruct, ier_flag

      IF (nextcur .eq. 0) nextcur = COUNT(extcur .ne. zero)
      WRITE (nwout, *) imse2 - 1, itse, nbsets, nobd, nextcur,
     1   nstore_seq
      IF (nbsets .gt. 0) WRITE (nwout, *) (nbfld(i),i=1,nbsets)
      WRITE (nwout, '(a)') mgrid_file

!-----------------------------------------------
!     Modification to obtain old fort.8 file
!-----------------------------------------------
      IF (loldout) THEN
        WRITE (nfort8, 710) gamma, nfp, ns, mpol, ntor, mnmax, itfsq,
     1         niter/100+1
!       the variables helsym (4th) and kpres (last) were set to 0.
        WRITE (nfort18, 40)voli,gamma,1.0/nfp,0.,mnmax,ns,mpol,ntor,
     1              ntor+1,1,itfsq,niter/100+1,0
      END IF
 40   FORMAT(1x,1p,e22.12,1x,3e12.5,8i4,i1)
 710  FORMAT(f10.3,7i6)

      IF (.not. lwrite) GOTO 1000

      pres(1) = pres(2)
      ntotal = 0

      DO js = 1, ns
         lj = (js - 1)*mnmax

         CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, js)
         write(23,*) js,(rmnc(mn),mn=1,mnmax)

         mn = 0
         bmod = SQRT(two*ABS(bsq(js,:nznt)-pres(js)))
         DO m = 0, mpol1
            nmin0 = -ntor
            IF (m .eq. 0) nmin0 = 0
            DO n = nmin0, ntor
               mn = mn + 1
               ntotal = MAX(mn,ntotal)
               dmult = two/(mscale(m)*nscale(ABS(n)))
               IF (m .eq. 0 .and. n .eq. 0) dmult = p5*dmult
               n1 = ABS(n)
               isgn = SIGN(1, n)
               gmn = 0
               bmn = 0
               currvmn = 0
               bsubumn = 0
               bsubvmn = 0
               bsubsmn = 0
               bsupumn = 0
               bsupvmn = 0
               lk = 0
               DO j = 1, ntheta2
                  DO k = 1, nzeta
                     lk = lk + 1
                     tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                         isgn*sinmui(j,m)*sinnv(k,n1))
                     tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                         isgn*cosmui(j,m)*sinnv(k,n1))
                     bmn = bmn + tcosi*bmod(lk)
                     gmn = gmn + tcosi*gsqrt(js,lk)
                     currvmn = currvmn + tcosi*izeta(js,lk)
                     bsubumn = bsubumn + tcosi*bsubu(js,lk)
                     bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
                     bsubsmn = bsubsmn + tsini*bsubs(js,lk)
                     bsupumn = bsupumn + tcosi*bsupu(js,lk)
                     bsupvmn = bsupvmn + tcosi*bsupv(js,lk)
                  END DO
               END DO
               IF (lasym) THEN
                  gmns = 0
                  bmns = 0
                  bsubumns = 0
                  bsubvmns = 0
                  bsubsmnc = 0
                  bsupumns = 0
                  bsupvmns = 0
                  lk = 0
                  DO j = 1, ntheta2
                     DO k = 1, nzeta
                        lk = lk + 1
                        tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                            isgn*sinmui(j,m)*sinnv(k,n1))
                        tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                            isgn*cosmui(j,m)*sinnv(k,n1))
                        bmns = bmns + tsini*bmod(lk)
                        gmns = gmns + tsini*gsqrt(js,lk)
                        bsubumns = bsubumns + tsini*bsubu(js,lk)
                        bsubvmns = bsubvmns + tsini*bsubv(js,lk)
                        bsubsmnc = bsubsmnc + tcosi*bsubs(js,lk)
                        bsupumns = bsupumns + tsini*bsupu(js,lk)
                        bsupvmns = bsupvmns + tsini*bsupv(js,lk)
                     END DO
                  END DO
               END IF
               IF (js .eq. ns/2) bmodmn(mn) = bmn
               IF (js .eq. ns) bmodmn1(mn) = bmn
               IF(js .eq. 1) THEN
                  IF (m .ne. NINT(xm(mn))) STOP 'm != NINT(xm)'
                  IF (n*nfp .ne. NINT(xn(mn))) STOP ' n = NINT(xn)'
                  WRITE (nwout, *) m, n*nfp
                  raxis_cc(0:ntor) = rmnc(1:ntor+1)
                  zaxis_cs(0:ntor) = zmns(1:ntor+1)
                  raxis_cs(0:ntor) = rmns(1:ntor+1)
                  zaxis_cc(0:ntor) = zmnc(1:ntor+1)
!                 ZERO HALF-MESH QUANTITIES
                  bsubumn = 0; bsubumns = 0
                  bsubvmn = 0; bsubvmns = 0
                  bsupumn = 0; bsupumns = 0
                  bsupvmn = 0; bsupvmns = 0
                  bmn = 0;     bmns     = 0
                  gmn = 0;     gmns     = 0
               END IF
!
!       PUT LMNS, LMNC ONTO HALF-MESH FOR BACKWARDS CONSISTENCY
!
               IF (js .eq. 1) THEN
                  lmnsh = 0
                  lmnch = 0
               ELSE IF (MOD(m,2) .eq. 0) THEN
                  lmnsh = p5*(lmns(mn) + lmnsjm1(mn))
                  lmnch = p5*(lmnc(mn) + lmncjm1(mn))
               ELSE
                  IF (js.eq.2 .and. m.eq.1) THEN
                     lmnsjm1(mn) = lmns(mn)
                     lmncjm1(mn) = lmnc(mn)
                  END IF
                  lmnsh = p5*(sm(js)*lmns(mn) + sp(js-1)*lmnsjm1(mn))
                  lmnch = p5*(sm(js)*lmnc(mn) + sp(js-1)*lmncjm1(mn))
               END IF

               WRITE (nwout, *) rmnc(mn), zmns(mn), lmnsh,
     1            bmn, gmn, bsubumn,
     2            bsubvmn, bsubsmn, bsupumn, bsupvmn, currvmn
               IF (lasym) THEN
                  WRITE (nwout, *) rmns(mn), zmnc(mn), lmnch,
     1            bmns, gmns, bsubumns, bsubvmns, bsubsmnc, 
     2            bsupumns, bsupvmns
               ENDIF

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
               IF (loldout) THEN
                 IF (js .eq. 1)
     1              WRITE (nfort8, 721) NINT(xm(mn)),NINT(xn(mn))
                 WRITE (nfort18,50)  xm(mn),xn(mn),rmnc(mn),zmns(mn),gmn
                 WRITE (nfort8, 731) rmnc(mn), zmns(mn), lmns(mn), bmn,
     1            gmn, bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
               END IF
            END DO
         END DO
         lmnsjm1 = lmns        !!Store previous full point values for averaging
         lmncjm1 = lmnc
      END DO
 50   FORMAT(1x,1p,5e14.6)
 721  FORMAT(2i10)
 731  FORMAT(5e20.13)

!
!     HALF-MESH QUANTITIES (except phi, jcuru, jcurv which are FULL MESH - computed in eqfor)
!     NOTE: jcuru, jcurv are local current densities, NOT integrated in s and normed to twopi
!     NOTE: In version <= 6.00, mass, press are written out in INTERNAL units
!     and should be multiplied by 1/mu0 to transform to pascals. In version > 6.00,
!     the pressure, mass are in correct (physical) units
!

      WRITE (nwout, *) (iotaf(js), presf(js)/mu0, phipf(js), phi(js), 
     1   jcuru(js)/mu0, jcurv(js)/mu0, js=1,ns)
      WRITE (nwout, *) (iotas(js), mass(js)/mu0, pres(js)/mu0,
     1   beta_vol(js), phip(js), buco(js), bvco(js), vp(js),
     2   overr(js), specw(js),js=2,ns)
!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      IF (loldout) THEN
         WRITE (nfort8, 731) (iotas(js),mass(js),pres(js),phips(js),
     1          buco(js),bvco(js),phi(js),vp(js),jcuru(js)/mu0,
     2          jcurv(js)/mu0,specw(js),js=2,ns)
        WRITE (nfort18,50)(iotas(js),mass(js),pres(js),-phips(js),
     1          vp(js),js=1,ns)
      END IF
!-----------------------------------------------

      WRITE (nwout, *) aspect, betatot, betapol, betator, betaxis, b0

!-----------------------------------------------
!     New output added to version 6.20
!-----------------------------------------------
      WRITE (nwout, *) NINT(signgs)
      WRITE (nwout, '(a)') input_extension
      WRITE (nwout, *) IonLarmor, VolAvgB, rbtor0, rbtor, ctor/mu0,
     1  Aminor_p, Rmajor_p, volume_p
!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      WRITE (nwout, *) (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     1       Dgeod(js), equif(js), js=2,ns-1)

      IF (lfreeb .and. nextcur .gt. 0) THEN
         WRITE (nwout, *) (extcur(i),i=1,nextcur)
         WRITE (nwout, *) (curlabel(i),i=1,nextcur)
      ENDIF

!-----------------------------------------------
!     NOTE: jdotb is in units of A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
!-----------------------------------------------
      WRITE (nwout, *) (fsqt(i),wdot(i),i=1,nstore_seq)
      WRITE (nwout, *) (jdotb(js),bdotgradv(js),js=1,ns)

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      IF (loldout) WRITE (nfort8, 731) (fsqt(i),wdot(i),i=1,100)

!-----------------------------------------------
!   Write diagno file
!-----------------------------------------------
      IF (lasym) THEN
!        DO mn = 1, mnmax
            potsin = 0; potcos = 0
!           DO mn0 = 1, mnpd
!              IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.
!    1              (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
!                 potsin = potvac(mn0)
!                 potcos = potvac(mn0+mnpd)
!                 EXIT
!              END IF
!           END DO
!        END DO

      ELSE
!        DO mn = 1, mnmax
            potsin = 0
!           DO mn0 = 1, mnpd
!              IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.
!    1              (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
!                 potsin = potvac(mn0)
!                 EXIT
!              END IF
!           END DO
!        END DO
         IF(ldiagno)THEN
           IF(lfreeb)THEN
             write(6,*)"Writing old diagno input file!"
             open(unit=21,file='diagno_input.data',status='unknown',
     1            action='write')
             write(21,'(a8)') "vmec2000"
             write(21,*) "nfp  mpol  ntor"
             write(21,*) nfp, mpol, ntor
             CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, ns)
             write(21,*) "rmnc"
             write(21,'(1p,e24.16)') (rmnc(mn),mn=1,mnmax)
             write(21,*) "zmns"
             write(21,'(1p,e24.16)') (zmns(mn),mn=1,mnmax)
             write(21,*) "potsin"
             DO mn0 = 1, mnpd
               write(21,'(1p,e24.16)') potvac(mn0)
             END DO
             write(21,*) "phiedge"
             write(21,*) phiedge
             write(21,*) "nextcur"
             write(21,*) nextcur
             write(21,*) "external currents"
             write(21,*) extcur(1:nextcur)
             CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc,1)
             write(21,*) "plasma current"
             write(21,*) ctor
             write(21,*) "plasma current filament fc R"
             write(21,*) rmnc(1:ntor+1)
             write(21,*) "plasma current filament fc z"
             write(21,*) zmns(1:ntor+1)
             close(unit=21)
           ELSE
             write(6,*)"Diagno-file request not completed!"
             write(6,*)"VMEC2000 not running in free-boundary mode!"
             write(6,*)"Check mgrid-file and flags!"
           ENDIF
         ENDIF
      END IF

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      IF (.not.lrecon) GOTO 900

      IF (imse2 - 1.gt.2 .or. itse.gt.0) THEN
         WRITE (nwout, *) tswgt, msewgt
         CALL smoothdata(nwout)

!       These knot values are on SQRT(s) grid
         presfactor = mu0*pthommax             !!*pfac moved to getthom
         WRITE (nwout, *) isnodes, (sknots(i),ystark(i),y2stark(i),i=
     1      1,isnodes)
         WRITE (nwout, *) ipnodes, (pknots(i),presfactor*ythom(i),
     1      presfactor*y2thom(i),i=1,ipnodes)
         WRITE (nwout, *)(datamse(i),rmid(i),qmid(i),shear(i),
     1      presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
         WRITE (nwout, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
         WRITE (nwout, *)(rthom(i),datathom(i),i=1,itse)
      ENDIF
      IF (nobd .gt. 0) THEN
         WRITE (nwout, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
         WRITE (nwout, *) flmwgt
      ENDIF
      IF (nbfldn .gt. 0) THEN
         DO n = 1, nbsets
            WRITE (nwout, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1         i=1,nbfld(n))
         END DO
         WRITE (nwout, *) bcwgt
      ENDIF

      WRITE (nwout, *) phidiam, delphid
!
!     Write Limiter & Prout plotting specs
!
      WRITE (nwout, *) nsets, nparts, nlim
      WRITE (nwout, *) (nsetsn(i),i=1,nsets)
      WRITE (nwout, *) (((pfcspec(i,j,k),i=1,nparts),j=1,nsetsn(k)),
     1   k=1,nsets)
      WRITE (nwout, *) (limitr(i), i=1,nlim)
      WRITE (nwout, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      WRITE (nwout, *) nrgrid, nzgrid
      WRITE (nwout, *) tokid
      WRITE (nwout, *) rx1, rx2, zy1, zy2, condif
      WRITE (nwout, *) imatch_phiedge
      IF (imatch_phiedge .eq. 2) CALL wroutlim(nwout)

 900  CONTINUE
!-----------------------------------------------
!     FREE BOUNDARY DATA
!-----------------------------------------------
!DEC$ IF DEFINED (NETCDF)
!DEC$ ELSE
      CALL freeb_data(rmnc, zmns, rmns, zmnc, bmodmn, bmodmn1)
!DEC$ ENDIF



 1000 CONTINUE

      WRITE (nwout, *) mgrid_mode

      CLOSE (nwout)

      END SUBROUTINE wrout
