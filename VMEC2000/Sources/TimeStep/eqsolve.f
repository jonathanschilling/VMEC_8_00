      SUBROUTINE eqsolve(nsval, interp_flag, ier_flag,
     1    lreset, lfirst, lscreen, reset_file_name)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, ns4, bad_init_jac_flag,
     1                       bad_jacobian_flag
      USE realspace
      USE vsvd
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nsval, ier_flag
      CHARACTER*(*) :: reset_file_name
      LOGICAL :: interp_flag, lreset, lfirst, lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p98 = 0.98_dp, p96 = 0.96_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsold=0, neqs2_old=0
      REAL(rprec) :: w1, r00s, w0, res0, wdota, r0dot
      REAL(rprec) :: delt0
      LOGICAL :: liter_flag, lreset_internal
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        hs      radial mesh size increment
!        iequi   counter used to call -EQFOR- at end of run
!        irzloff offset in xc array between R,Z,L components
!        ijacob  counter for number of times jacobian changes sign
!        irst    counter monitoring sign of jacobian; resets R, Z, and
!                Lambda when jacobian changes sign and decreases time step
!        signgs  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
!        iterj   stores position in main iteration loop (j=1,2)
!        itfsq   counter for storing FSQ into FSQT for plotting
!        ivac    counts number of free-boundary iterations
!        ndamp   number of iterations over which damping is averaged
!        meven   parity selection label for even poloidal modes of R and Z
!        modd    parity selection label for odd poloidal modes of R and
!        gc      stacked array of R, Z, Lambda Spectral force coefficients (see readin for stack order)
!        xc      stacked array of scaled R, Z, Lambda Fourier coefficients

!
!       INITIALIZE MESH-DEPENDENT SCALARS
!
      ns = nsval
      ns1 = ns - 1
      hs = one/ns1
      ohs = one/hs
      mns = ns*mnsize
      irzloff = ntmax*mns
      nrzt = nznt*ns
      neqs = 3*irzloff
      neqs1 = neqs + 1
      neqs2 = neqs1 + 1
      ijacob = 0
      delt0   = delt

      WRITE (nthreed, 10) ns, mnmax, ftolv
      IF (lscreen) PRINT 10, ns, mnmax, ftolv
   10 FORMAT(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',
     1   1p,e10.3)
      IF (imovephi .gt. 0 .and. lscreen) PRINT *,
     1   'Turning on Ip Matching by varying Phi-Edge.'

!
!     ALLOCATE NS-DEPENDENT ARRAYS
!
      CALL allocate_ns(lreset, interp_flag, neqs2_old)

      lreset_internal = (ns .gt. nsold)

      IF (lreset_internal .and. neqs2_old.gt.0) gc(1:neqs2_old) =
     1  scalxc(1:neqs2_old)*xstore(1:neqs2_old)

      lreset_internal = lreset_internal .or. lreset

 1000 CONTINUE

      itfsq = 0
      fsq     = one
      rsfac   = one
      w1      = zero
      r00s    = zero
      gphifac = zero
      grmse   = zero

!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
   20 CONTINUE
      iter2 = 1
      irst = 1
      CALL profil1d (xc, xcdot, lreset)
      IF (lfirst .and. .not.lreset)               !!Command line with lreset=F
     1    CALL load_xc_from_wout(xc(1), xc(1+irzloff), xc(1+2*irzloff),
     2         ntor, mpol1, ns, reset_file_name, lreset_internal)
      CALL profil3d (xc(1), xc(1+irzloff), lreset_internal)

!
!     INTERPOLATE FROM COARSE TO NEXT FINER RADIAL GRID
!
      IF (interp_flag) CALL interp (xc, gc, scalxc, ns, nsold)
      nsold = ns
      neqs2_old = neqs2

!
!     Store XC,XCDOT for possible restart
!
      CALL restart_iter(delt0)
      iter1 = iter2
      liter_flag = .true.
      ier_flag = 0

!
!     ENTER FORCE ITERATION LOOP
!
      iter_loop: DO WHILE (liter_flag)
!
!     ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
!
         CALL evolve (delt0, ier_flag, liter_flag, lscreen, r0dot)
         IF (ijacob.eq.0 .and. ier_flag.eq.bad_init_jac_flag
     1                   .and. ns.ge.3) THEN
            IF (lscreen) PRINT *,
     1         ' TRYING TO IMPROVE INITIAL MAGNETIC AXIS GUESS'
            CALL guess_axis (r1, z1, ru0, zu0)
            lreset_internal = .true.
            ijacob = 1
            GOTO 20
         ELSE IF (ier_flag .ne. 0) THEN
            IF (ier_flag .eq. bad_init_jac_flag) 
     1          ier_flag = bad_jacobian_flag
            RETURN
         ENDIF

         w0 = wb + wp/(gamma - one)

!
!     ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)
!
         IF (ijacob .eq. 25) THEN
            irst = 2
            CALL restart_iter(delt0)
            delt0 = p98*delt
            IF (lscreen) PRINT 120, delt0
            GOTO 1000
         ELSE IF (ijacob .eq. 50) THEN
            irst = 2
            CALL restart_iter(delt0)
            delt0 = p96*delt
            IF (lscreen) PRINT 120, delt0
            GOTO 1000
         ELSE IF (ijacob .ge. 75) THEN
            ier_flag = 2
            liter_flag = .false.
         ELSE IF (iter2.ge.niter .and. liter_flag) THEN
            ier_flag = 4
            liter_flag = .false.
         ENDIF

!
!       TIME STEP CONTROL
!
         IF (iter2 .eq. iter1) res0 = fsq
         res0 = MIN(res0,fsq)
!       Store current state (irst=1)
         IF (fsq .le. res0 .and. iter2-iter1 .gt. 10) THEN
            CALL restart_iter(delt0)
!       Residuals are growing in time, reduce time step
         ELSE IF (fsq .gt. 100.0_dp*res0 .and. iter2 .gt. iter1) THEN
            irst = 2
         ELSE IF (iter2 - iter1 .gt. ns4/2 .and. iter2 .gt. 2*ns4
     1        .and. fsqr+fsqz .gt. c1pm2) THEN
            irst = 3
         ENDIF

         IF (irst .ne. 1) THEN
!       Retrieve previous good state
            CALL restart_iter(delt0)
            iter1 = iter2
         ELSE
!       Increment time step and Printout every nstep iterations
            IF (MOD(iter2,nstep).eq.0 .or. iter2.eq.1 .or.
     1         .not.liter_flag) CALL printout(iter2, delt0, w0, lscreen)
            iter2 = iter2 + 1
         ENDIF

!       Store force residual, wdot for plotting
         wdota = ABS(w0 - w1)/w0
         r0dot = ABS(r00 - r00s)/r00
         r00s = r00
         w1 = w0
         IF (ivac.eq.1 .and. lreset) THEN
            IF (lscreen) PRINT 110, iter2
            WRITE (nthreed, 110) iter2
            ivac = ivac + 1
         ENDIF
!
!       STORE FSQ FOR PLOTTING. EVENTUALLY, STORE FOR EACH RADIAL MESH
!
         IF (MOD(iter2,niter/nstore_seq + 1).eq.0 .and. ns.eq.
     1      ns_array(multi_ns_grid)) THEN
            IF (itfsq .lt. nstore_seq) THEN
              itfsq = itfsq + 1
              fsqt(itfsq) = fsqr + fsqz
              wdot(itfsq) = MAX(wdota,c1pm13)
            END IF
         END IF

      END DO iter_loop

      WRITE (nthreed, 60) wdota, r0dot, w0*twopi**2
      IF (lrecon) WRITE (nthreed, 70) r00*fsqsum0/wb

   60 FORMAT(/,' d(ln W)/dt = ',1p,e14.3,' d(ln R0)/dt = ',e14.3,/,
     1         ' MHD Energy = ',e14.6)
   70 FORMAT(' Average radial force balance: Int[FR(m=0)]',
     1   '/Int(B**2/R) = ',1p,e12.5,' (should tend to zero)'/)
  110 FORMAT(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)
  120 FORMAT(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,
     1  /,2x,'If this does NOT resolve problem, try changing ',
     2       '(decrease OR increase) the value of DELT')

      END SUBROUTINE eqsolve
