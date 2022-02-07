      SUBROUTINE readin(input_file, iseq_count, lfirst, ier_flag,
     1      lscreen)
      USE vmec_main
      USE vmec_params
      USE vacmod
      USE vsvd
      USE vspline
      USE timer_sub
      USE mgrid_mod, ONLY: nextcur, curlabel, nfper0, read_mgrid
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER iseq_count, ier_flag
      LOGICAL lfirst, lscreen
      CHARACTER*(*) :: input_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: ns_default = 31
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iexit, ipoint, n, iunit, ier_flag_init,
     1   i, ni, m, nsmin, igrid, mj, isgn, ioff, joff
      REAL(rprec), DIMENSION(:,:), POINTER ::
     1  rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      REAL(rprec) :: rtest, ztest, tzc, trc, delta
      CHARACTER*(100) :: line, line2
C-----------------------------------------------
!
!       LOCAL VARIABLES
!
!       rbcc,rbss,rbcs,rbsc
!                boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
!       zbcc,zbss,zbcs,zbsc
!                boundary Fourier coefficient arrays for Z
!
!       XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!
!       STACKING ORDER DEPENDS ON LASYM AND LTHREED. EACH COMPONENT XCC, XSS, XSC, XCS
!       HAS SIZE = mns. (PHIFAC, MSE TAKE UP 1 INDEX EACH AT END OF ARRAY)
!
!         LTHREED=F,      LTHREED=F,      LTHREED=T,      LTHREED=T 
!         LASYM=F         LASYM=T         LASYM=F         LASYM=T
!
!          rmncc           rmncc           rmncc           rmncc           
!          zmnsc           rmnsc           rmnss           rmnss
!          lmnsc           zmnsc           zmnsc           rmnsc
!                          zmncc           zmncs           rmncs
!                          lmnsc           lmnsc           zmnsc
!                          lmncc           lmncs           zmncs
!                                                          zmncc
!                                                          zmnss
!                                                          lmnsc
!                                                          lmncs
!                                                          lmncc
!                                                          lmnss
!
!
!                STANDARD INPUT DATA AND RECOMMENDED VALUES
!
!   Plasma parameters (MKS units)
!          ai:   iota (ncurr=0) or toroidal current density (=ac, ncurr=1)
!                expansion coefficients (series in s)
!          am:   mass or pressure (gamma=0) expansion coefficients (series in s)
!                in MKS units (NWT/M**2)
!      curtor:   value of toroidal current (A). Used if ncurr = 1 to specify
!                current profile, or IF in data reconstruction mode.
!      extcur:   array of currents in each external current group. Used to
!                multiply Green''s function for fields and loops read in from
!                MGRID file. Should use real current units (A).
!       gamma:   value of compressibility index (gamma=0 => pressure prescribed)
!         nfp:   number of toroidal field periods ( =1 for Tokamak)
!         rbc:   boundary coefficients of COS(m*theta-n*zeta) for R
!         zbs:   boundary coefficients of SIN(m*theta-n*zeta) for Z
!         rbs:   boundary coefficients of SIN(m*theta-n*zeta) for R
!         zbc:   boundary coefficients of COS(m*theta-n*zeta) for Z
!
!
!   Numerical parameters
!  mgrid_file:   full path for vacuum green''s FUNCTION data
!       ncurr:   flux conserving (=0) or prescribed toroidal current (=1)
!    ns_array:   array of radial mesh SIZEs to be used in multigrid sequence
!
!   Convergence control parameters
!  ftol_array:   array of value of residual(s) at which each multigrid
!                iteration ENDs
!       niter:   number of iterations (used to terminate run)
!       nstep:   number of timesteps between printouts on screen
!    nvacskip:   iterations skipped between full update of vacuum solution
!       tcon0:   weight factor for constraint force (=1 by DEFAULT)
!
!   Equilibrium reconstruction parameters
!      phifac:   factor scaling toroidal flux to match apres or limiter
!   datastark:   pitch angle data from stark measurement
!    datathom:   pressure data from Thompson, CHEERS (Pa)
!     imatch_         = 1 (default),match value of PHIEDGE in input file
!     phiedge:   = 0, USE pressure profile width to determine PHIEDGE
!                = 2, USE LIMPOS data (in mgrid file) to find PHIEDGE
!                = 3, USE Ip to find PHIEDGE (fixed-boundary only)
!        imse:   number of Motional Stark effect data points
!                >0, USE mse data to find iota; <=0, fixed iota profile ai
!        itse:   number of pressure profile data points
!                = 0, no thompson scattering data to READ
!     isnodes:   number of iota spline points (computed internally unless specified explicitly)
!     ipnodes:   number of pressure spline points (computed internally unless specified explicitly)
!       lpofr:   LOGICAL variable. =.true. IF pressure data are
!                prescribed in REAL space. =.false. IF data in flux space.
!      pknots:   array of pressure knot values in SQRT(s) space
!      sknots:   array of iota knot values in SQRT(s) space
!       tensp:   spline tension for pressure profile
!
!       tensi:   spline tension for iota
!      tensi2:   vbl spline tension for iota
!      fpolyi:   vbl spline tension form factor (note: IF tensi!=tensi2
!               THEN tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
!               - - - - - - - - - - - - - - - - - -
!    mseangle_   uniform EXPerimental offset of MSE data
!     offset:    (calibration offset) ... PLUS ...
!    mseangle_   multiplier on mseprof offset array
!     offsetM:   (calibration offset)
!     mseprof:   offset array from NAMELIST MSEPROFIL
!                so that the total offset on the i-th MSE data point is
!                taken to be
!                = mseangle_offset+mseangle_offsetM*mseprof(i)
!               - - - - - - - - - - - - - - - - - -
! pres_offset:   uniform arbitrary  radial offset of pressure data
!     presfac:   number by which Thomson scattering data is scaled
!                to get actual pressure
!     phidiam:   diamagnetic toroidal flux (Wb)
!      dsiobt:   measured flux loop signals corresponding to the
!                combination of signals in iconnect array
!     indxflx:   array giving INDEX of flux measurement in iconnect array
!    indxbfld:   array giving INDEX of bfield measurement used in matching
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!    bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
!                the orientation br*COS(abcoil) + bz*SIN(abcoil)
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!                of the m-th coil in the n-th set from mgrid file.
!       nflxs:   number of flux loop measurements used in matching
!    nbfld(n):   number of selected EXTERNAL bfield measurements in set n from nml file
!      nbfldn:   total number of EXTERNAL bfield measurements used in matching
!               - - - - - - - - - - - - - - - - - -
!             NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
!             AS PERCENT OF RESPECTIVE MEASUREMENT
!  sigma_thom:   standard deviation (Pa) for pressure profile data
! sigma_stark:   standard deviation (degrees) in MSE data
!  sigma_flux:   standard deviaton (Wb) for EXTERNAL poloidal flux data
!     sigma_b:   standard deviation (T) for EXTERNAL magnetic field data
!sigma_current:  standard deviation (A) in toroidal current
!sigma_delphid:  standard deviation (Wb) for diamagnetic match
!
!
!       THE (ABSOLUTE) CHI-SQ ERROR IS DEFINED AS FOLLOWS:
!
!          2
!       CHI      =     SUM [ EQ(K,IOTA,PRESSURE)  -  DATA(K) ] ** 2
!                     (K) -----------------------------------
!                                   SIGMA(K)**2
!
!       HERE, SIGMA IS THE STANDARD DEVIATION OF THE MEASURED DATA, AND
!       EQ(IOTA,PRESSURE) IS THE EQUILIBRIUM EXPRESSION FOR THE DATA TO BE
!       MATCHED:
!
!       EQ(I)   =    SUM [ W(I,J)*X(J) ]
!                   (J)
!
!       WHERE W(I,J) ARE THE (LINEAR) MATRIX ELEMENTS AND X(J) REPRESENT
!       THE KNOT VALUES OF IOTA (AND/OR PRESSURE). THE RESULTING LEAST-SQUARES
!       MATRIX ELEMENTS AND DATA ARRAY CAN BE EXPRESSED AS FOLLOWS:
!
!       ALSQ(I,J) = SUM [ W(K,I) * W(K,J) / SIGMA(K) ** 2]
!                   (K)
!
!       BLSQ(I)   = SUM [ W(K,I) * DATA(K)/ SIGMA(K) ** 2]
!                   (K)
!
!       THEREFORE, INTERNALLY IT IS CONVENIENT TO WORK WITH THE 'SCALED'
!       W'(K,I) = W(K,I)/SIGMA(K) AND DATA'(K) = DATA(K)/SIGMA(K)
!
!       ****!   I - M - P - O - R - T - A - N - T     N - O - T - E   *****
!
!       THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
!       SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
!       TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
!       AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
!       THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
!       SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.
!
      ier_flag_init = ier_flag
      ier_flag = 0
      CALL second0(treadon)
      IF (ier_flag_init .eq. 4) GOTO 1000

!
!     READ IN DATA FROM INDATA FILE
!
      CALL read_indata(input_file, iunit, ier_flag)
      IF (ier_flag .ne. 0) RETURN


      IF (tensi2 .eq. zero ) tensi2 = tensi

!
!     Open output files here, print out heading to threed1 file
!
      CALL heading(input_extension, time_slice,
     1      iseq_count, lmac, lscreen, lfirst)

!
!     READ IN COMMENTS DEMARKED BY "!"
!
      REWIND (iunit)
      iexit = 0
      DO WHILE( iexit.eq.0 )
         READ (iunit, '(a)') line
         iexit = INDEX(line,'INDATA') + index(line,'indata')
         ipoint = INDEX(line,'!')
         IF (ipoint .eq. 1) WRITE (nthreed, *) TRIM(line)
      ENDDO
      CLOSE (iunit)

!
!     READ IN AND STORE (FOR SEQUENTIAL RUNNING) MAGNETIC FIELD DATA
!     FROM MGRID_FILE FIRST TIME (lfirst = T) ONLY
!     SET LOGICAL FLAGS FOR ALL SUBSEQUENT RUNS
!
      IF (lfirst .and. lfreeb) THEN
         CALL second0(trc)
         CALL read_mgrid (mgrid_file, extcur, nv, nfp, 
     1                    lscreen, ier_flag)
         CALL second0(tzc)

         IF (lfreeb) THEN
            WRITE (6,'(2x,a,1p,e10.2,a)') 'Time to read MGRID file: ', 
     1             tzc - trc, ' s'
            IF (ier_flag .ne. 0) RETURN
            WRITE (nthreed,20) nr0b, nz0b, np0b, rminb, rmaxb,
     1                         zminb, zmaxb, TRIM(mgrid_file)
 20         FORMAT(//,' VACUUM FIELD PARAMETERS:',/,1x,24('-'),/,
     1     '  nr-grid  nz-grid  np-grid      rmin      rmax      zmin',
     2     '      zmax     input-file',/,3i9,4f10.3,5x,a)
         END IF
      END IF

!
!     PARSE NS_ARRAY
!
      nsin = MAX (3, nsin)
      multi_ns_grid = 1
      IF (ns_array(1) .eq. 0) THEN                   !!Old input style
          ns_array(1) = MIN(nsin,nsd)
          multi_ns_grid = 2
          ns_array(multi_ns_grid) = ns_default        !!Run on 31-point mesh
      ELSE
          nsmin = 1
          DO WHILE (ns_array(multi_ns_grid) .gt. nsmin .and.
     1             multi_ns_grid .lt. 100)
             nsmin = MAX(nsmin, ns_array(multi_ns_grid))
            IF (nsmin.le.nsd) THEN
               multi_ns_grid = multi_ns_grid + 1
            ELSE                                      !!Optimizer, Boozer code overflows otherwise
               ns_array(multi_ns_grid) = nsd
               nsmin = nsd
               PRINT *,' NS_ARRAY ELEMENTS CANNOT EXCEED ',nsd
               PRINT *,' CHANGING NS_ARRAY(',multi_ns_grid,') to ', nsd
            END IF
          END DO
          multi_ns_grid = multi_ns_grid - 1
      ENDIF
      IF (ftol_array(1) .eq. zero) THEN
         ftol_array(1) = 1.e-8_dp
         IF (multi_ns_grid .eq. 1) ftol_array(1) = ftol
         DO igrid = 2, multi_ns_grid
            ftol_array(igrid) = 1.e-8_dp * (1.e8_dp * ftol)**
     1        ( REAL(igrid-1,rprec)/(multi_ns_grid-1) )
         END DO
      ENDIF

!
!     WRITE OUT DATA TO THREED1 FILE
!
      WRITE (nthreed,100)
     1  ns_array(multi_ns_grid),ntheta1,nzeta,mpol,ntor,nfp,
     2  gamma,spres_ped,phiedge,curtor
 100  FORMAT(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'/,
     1  1x,45('-'),/,
     2  '     ns     nu     nv     mu     mv',/,
     3  5i7,//,' CONFIGURATION PARAMETERS:',/,1x,25('-'),/,
     4  '    nfp      gamma      spres_ped    phiedge(wb)     curtor(A)'
     5  ,/,i7,1p,e11.3,2e15.3,e14.3,/)

      IF (nvacskip.le.0) nvacskip = nfp
      WRITE (nthreed,110) ncurr,niter,ns_array(1),nstep,nvacskip,
     1  ftol_array(multi_ns_grid),tcon0,lasym
 110    FORMAT(' RUN CONTROL PARAMETERS:',/,1x,23('-'),/,
     1  '  ncurr  niter   nsin  nstep  nvacskip      ftol     tcon0',
     2  '   lasym',/, 4i7,i10,1p,2e10.2,L8/)

      IF (nextcur .gt. 0) THEN
         WRITE(nthreed, "(' EXTERNAL CURRENTS',/,1x,17('-'))")
         ni = 0
         IF (ALLOCATED(curlabel))
     1      ni = MAXVAL(LEN_TRIM(curlabel(1:nextcur)))
         ni = MAX(ni+4, 14)
         WRITE (line,  '(a,i2.2,a)') "(5a",ni,")"
         WRITE (line2, '(a,i2.2,a)') "(5(",ni-12,"x,1p,e12.4))"
         DO i = 1,nextcur,5
            ni = MIN(i+4, nextcur)
            IF (ALLOCATED(curlabel))
     1      WRITE (nthreed, line, iostat=mj) (TRIM(curlabel(n)),n=i,ni)
            WRITE (nthreed, line2,iostat=mj) (extcur(n), n=i,ni)
         ENDDO
         WRITE (nthreed, *)
      ENDIF

      IF (bloat .ne. one) THEN
          WRITE (nthreed,'(" Profile Bloat Factor: ",f7.5)') bloat
          phiedge = phiedge*bloat
      ENDIF

      IF (pres_scale .ne. one) THEN
          WRITE (nthreed,'(" Pressure profile factor: ",f7.5,
     1           " (embedded in mass)")') pres_scale
          am = am * pres_scale
      END IF

      WRITE(nthreed,130)ipmass
 130  FORMAT(' MASS PROFILE COEFFICIENTS am - newton/m**2',
     1  ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'),'IPMASS=',i3)
      WRITE(nthreed,135)(am(i-1),i=1, SIZE(am))
      IF (ncurr.eq.0) THEN
          WRITE(nthreed,140)ipiota
          WRITE(nthreed,135)(ai(i-1),i=1, SIZE(ai))
      ELSE
          WRITE(nthreed,145)ipcurr
          WRITE(nthreed,135)(ac(i-1),i=1, SIZE(ac))
      ENDIF
!     WRITE(nthreed,150)
!     WRITE(nthreed,135)(aphi(i-1),i=1, SIZE(aphi))    !!Temporarily disabled

 135  FORMAT(1p,6e12.3)
 140  FORMAT(/' IOTA PROFILE COEFFICIENTS ai',
     1   ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'),'IPIOTA=',i3)
 145  FORMAT(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1 ' ac (EXPANSION IN TOROIDAL FLUX):',/,1x,38('-'),'IPCURR=',i3)
 150  FORMAT(/' NORMALIZED TOROIDAL FLUX COEFFICIENTS aphi',
     1   ' (EXPANSION IN S):',/,1x,35('-'))
      WRITE(nthreed,180)
 180  FORMAT(/,' R-Z FOURIER BOUNDARY COEFFICIENTS AND',
     1         ' MAGNETIC AXIS INITIAL GUESS',/,
     1  ' R = RBC*COS(m*u - n*v) + RBS*SIN(m*u - n*v),',
     2  ' Z = ZBC*COS(m*u - n*v) + ZBS*SIN(m*u-n*v)'/1x,86('-'),
     3  /,'   nb  mb     rbc         rbs         zbc         zbs   ',
     4   '    raxis(c)    raxis(s)    zaxis(c)    zaxis(s)')

 1000  CONTINUE

      IF (.not.lasym) THEN
!
!       CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
!

      delta = ATAN( (rbs(0,1) - zbc(0,1))/
     1           (ABS(rbc(0,1)) + ABS(zbs(0,1))) )
      IF (delta .ne. zero) THEN
        DO m = 0,mpol1
          DO n = -ntor,ntor
            trc = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
            rbs(n,m) = rbs(n,m)*COS(m*delta) - rbc(n,m)*SIN(m*delta)
            rbc(n,m) = trc
            tzc = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
            zbs(n,m) = zbs(n,m)*COS(m*delta) - zbc(n,m)*SIN(m*delta)
            zbc(n,m) = tzc
          ENDDO
        ENDDO
      ENDIF

      ENDIF

!
!     ALLOCATE MEMORY FOR NU, NV, MPOL, NTOR SIZED ARRAYS
!
      CALL allocate_nunv

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!         + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!         + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!
!
!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)
!     ARRAY STACKING ORDER DETERMINED HERE
!

      rbcc => rmn_bdy(:,:,rcc)
      zbsc => zmn_bdy(:,:,zsc)
      IF (lthreed) THEN
         rbss => rmn_bdy(:,:,rss)
         zbcs => zmn_bdy(:,:,zcs)
      END IF

      IF (lasym) THEN
         rbsc => rmn_bdy(:,:,rsc)
         zbcc => zmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => zmn_bdy(:,:,zss)
         END IF
      ENDIF

      rmn_bdy = 0;  zmn_bdy = 0

      ioff = LBOUND(rbcc,1)
      joff = LBOUND(rbcc,2)

      DO m=0,mpol1
         mj = m+joff
         DO n=-ntor,ntor
            ni = ABS(n) + ioff
            IF (n .eq. 0) THEN
               isgn = 0
            ELSE IF (n .gt. 0) THEN
               isgn = 1
            ELSE
               isgn = -1
            END IF
            rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
            zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)
            IF (m .eq. 0) zbsc(ni,mj) = 0

            IF (lthreed) THEN
               rbss(ni,mj) = rbss(ni,mj) + isgn*rbc(n,m)
               zbcs(ni,mj) = zbcs(ni,mj) - isgn*zbs(n,m)
               IF (m .eq. 0) rbss(ni,mj) = 0
            END IF

            IF (lasym) THEN
               rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
               zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)
               IF (m .eq. 0) rbsc(ni,mj) = 0
               IF (lthreed) THEN
               rbcs(ni,mj) = rbcs(ni,mj) - isgn*rbs(n,m)
               zbss(ni,mj) = zbss(ni,mj) + isgn*zbc(n,m)
               IF (m .eq. 0) zbss(ni,mj) = 0
               END IF
            END IF

            IF (ier_flag_init .ne. 0) CYCLE
            trc = ABS(rbc(n,m)) + ABS(rbs(n,m))
     1          + ABS(zbc(n,m)) + ABS(zbs(n,m))
            IF (m .eq. 0) THEN
               IF (n .lt. 0) CYCLE
               IF (trc.eq.zero .and. ABS(raxis_cc(n)).eq.zero .and.
     1             ABS(zaxis_cs(n)).eq.zero) CYCLE
               WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m), raxis_cc(n), raxis_cs(n),
     2                   zaxis_cc(n), zaxis_cs(n)
            ELSE
               IF (trc .eq. zero) CYCLE
               WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m)
            END IF
         END DO
      END DO
 195  FORMAT(i5,i4,1p,8e12.4)

!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!
      m = 1
      mj = m+joff
      rtest = SUM(rbcc(1:ntor1,mj))
      ztest = SUM(zbsc(1:ntor1,mj))
      signgs = one
      IF (rtest*ztest .gt. zero) signgs = -one


      iresidue = -1
      IF (lrecon) THEN
!
!       DETERMINE CURRENT-FLUX CONSISTENCY CHECK
!
        signiota = one
        IF (signgs*curtor*phiedge .lt. zero)signiota = -one
        IF (sigma_current .eq. zero) THEN
          WRITE (*,*) 'Sigma_current cannot be zero!'
          ier_flag = 1
          RETURN
        END IF

!
!       SET UP RECONSTRUCTION FIXED PROFILES
!
        dcon = ATAN(one)/45
        CALL readrecon                   !Setup for reconstruction mode
        CALL fixrecon(ier_flag)          !Fixed arrays for reconstruction
        IF (ier_flag .ne. 0) RETURN
      END IF

      currv = mu0*curtor              !Convert to Internal units

      CALL second0(treadoff)
      timer(tread) = timer(tread) + (treadoff-treadon)

      END SUBROUTINE readin
