c !DEC$ IF DEFINED (RISC)
c @PROCESS ALIAS(ARYOVRLP)
c !DEC$ ENDIF
      MODULE precon2d
      USE stel_kinds, ONLY: rprec2 => rprec, dp2 => dp
      USE vmec_dim
      USE vmec_params, mnyq_vp => mnyq, nnyq_vp => nnyq
      USE vparams, ONLY: nthreed, zero
      USE vmec_input, ONLY: ntor, nzeta, lfreeb, lasym, ncurr
      USE timer_sub
      USE fbal, ONLY: lforbal

      INTEGER, PARAMETER :: sp = KIND(1.0)
!      INTEGER, PARAMETER :: sp = rprec2
      INTEGER :: ntyptot, ntype_2d, m_2d, n_2d
      INTEGER, ALLOCATABLE :: ipiv_blk(:,:)
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) ::
     1    block_diag, block_plus, block_mins

      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) ::
     1    blockl_diag, blockl_plus, blockl_mins
      REAL(rprec2), DIMENSION(:,:,:,:), ALLOCATABLE :: gc_save
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)
     1             :: block_dsave, block_msave, block_psave
      REAL(rprec2), ALLOCATABLE :: rzl_save(:,:)
      REAL(rprec2), ALLOCATABLE, DIMENSION(:,:,:) ::
     1   r1_save, ru_save, rv_save, rcon_save,
     2   z1_save, zu_save, zv_save, zcon_save,
     3   lu_save, lv_save
      LOGICAL :: lprec2d = .false.
      LOGICAL :: loglam_blks=.true., lsweep_fast = .true.,
     1           ljog_test = .false., l_backslv = .false.
      PRIVATE :: swap_forces, reswap_forces, gc_save, block_dsave,
     1           block_msave, block_psave
!
!     SP:        forces single precision for blocks (smaller size)
!     LOGLAM_BLKS: internal variable used to control call to lamblks subroutine
!     LPREC2D:   used externally to control when block_precond is called
!                it is set (.true.) when compute_blocks is called
!     L_BACKSLV: if true, test that Hessian is inverted correctly by back-solving
!     LJOG_TEST: used to test analytic Hessian calculation against reliable
!                "jog" technique (when true). ONLY WORKS with lsweep_fast=.false.,
!                because lsweep_fast=true can not compute lambda forces (which are not symmetric)
!     LSWEEP_FAST : if true, uses sweep2 fast (symmetric R,Z) Hessian; if false (which is most
!                   general), uses sweep3, 2/3 slower Hessian
!     RZL_SAVE:  save offset of m=1 r,z modes at js=1 when iequi=2
!
      CONTAINS

      SUBROUTINE swap_forces(gc, temp, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(rprec2), DIMENSION(nblocks,mblk), INTENT(in)  :: gc
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(out) :: temp
C-----------------------------------------------
!
!     orders forces (gc) array prior to applying
!     block-tridiagonal pre-conditioner. on exit, temp is ordered
!     flip sign so eigenvalue is negative (corresponding to damping)
!
      temp = -TRANSPOSE(gc)

      END SUBROUTINE swap_forces

      SUBROUTINE reswap_forces(temp, gc, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(rprec2), DIMENSION(nblocks,mblk), INTENT(inout)  :: gc
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(in) :: temp
C-----------------------------------------------
!
!     Following application of block pre-conditioner, restores original
!     order of forces (gc) array previously ordered by call to "swap_forces"
!
      gc = TRANSPOSE(temp)

      END SUBROUTINE reswap_forces

      SUBROUTINE block_precond(gc)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mblk, m, n, js, ntype, istat
      REAL(rprec2), ALLOCATABLE, DIMENSION(:,:) :: temp
      REAL(rprec2) :: t1, error
C-----------------------------------------------
!
!     Applies 2D block-preconditioner to forces vector (gc)
!
      IF (ntyptot .le. 0) STOP 'ntyptot must be > 0'

      mblk = ntyptot*mnsize
      ALLOCATE (temp(mblk,ns), stat=ntype)
      IF (ntype .ne. 0) STOP 'Allocation error1 in block_precond'

      CALL swap_forces(gc, temp, mblk, ns)

!     Perform preconditioning, using LU factors stored in block_... matrices
      CALL blk3d_slv(block_diag, block_mins, block_plus, temp,
     1                  ipiv_blk, mblk, ns)

      CALL reswap_forces(temp, gc, mblk, ns)

      IF (l_backslv) THEN
         l_backslv = .false.

         WRITE (6, *) ' Writing block Hessian check to unit 34'
         WRITE (34, *)
         WRITE (34, *) ' BLK3D FACTORIZATION CHECK: Ax = b ?'

         DO n = 0, ntor
            WRITE (34, *) ' N = ', n
            DO m = 0, mpol1
               WRITE (34, *) ' M = ', m
               DO ntype = 1, ntyptot
                  WRITE (34, *) ' TYPE = ', ntype
                  js = 1
                  t1 = SUM(block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:)
     1               + block_psave(n,m,ntype,:,:,:,js)*gc(js+1,:,:,:))

                  error = t1 + gc_save(js,n,m,ntype)
                  IF (abs(error) .gt. 1.E-12)
     1            WRITE (34, *) ' js = ', js,' Ax = ', t1,' b = ',
     2                 -gc_save(js,n,m,ntype), 'Error = ', error

               DO js = 2, ns-1
                  t1 = SUM
     1                  (block_msave(n,m,ntype,:,:,:,js)*gc(js-1,:,:,:)
     2                 + block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:)
     3                 + block_psave(n,m,ntype,:,:,:,js)*gc(js+1,:,:,:))
               error = t1 + gc_save(js,n,m,ntype)
               IF (abs(error) .gt. 1.E-12)
     1         WRITE (34, *) ' js = ', js,' Ax = ', t1,' b = ',
     2                 -gc_save(js,n,m,ntype), 'Error = ', error
               END DO

               js = ns
               t1 = SUM(block_msave(n,m,ntype,:,:,:,js)*gc(js-1,:,:,:)
     1            + block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:))
               error = t1 + gc_save(js,n,m,ntype)
               IF (abs(error) .gt. 1.E-12)
     1         WRITE (34, *) ' js = ', js,' Ax = ', t1,' b = ',
     2                 -gc_save(js,n,m,ntype), 'Error = ', error
            END DO
         END DO
         END DO

         DEALLOCATE(block_dsave, block_msave, block_psave, gc_save,
     1              stat=istat)

      END IF

      DEALLOCATE (temp, stat=istat)

      END SUBROUTINE block_precond

      SUBROUTINE compute_blocks (xc, xcdot, gc)
      USE realspace, ONLY: sqrts
      USE vmec_main, ONLY: liota, iter2, iequi, lthreed
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2),DIMENSION(ns,0:ntor,0:mpol1,3*ntmax) :: xc, gc, xcdot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, i, iequi_save, ntype, istat,
     1           lamtype, rztype, ncurr_save, mblk,
     2           nmax_jog, ibsize, iunit
      REAL(rprec2) :: time_on, time_off, diag_val(ns), diag_lam(ns),
     1               bsize, tprec2don, tprec2doff
      CHARACTER*100 :: label
      LOGICAL, PARAMETER :: lscreen = .false.
C-----------------------------------------------
      DO i = 1,2
         IF (i .eq. 1) iunit = 6
         IF (i .eq. 2) iunit = nthreed
         WRITE (iunit, '(/,2x,2a,i5,a)') 'Initializing 2D block',
     1         ' preconditioner at ', iter2,' iterations'
      END DO

      CALL second0(tprec2don)

      lprec2d = .true.
      ntyptot = SIZE(gc,4)
      rztype = 2*ntmax
      lamtype = rztype+1
      lprec2d = .true.
      mblk = ntyptot*mnsize

      bsize = mblk*mblk*3*ns*KIND(block_diag)
      IF (bsize .lt. 1.E6_dp2) THEN
          ibsize = bsize/1.E3_dp2
          label = " Kb"
      ELSE IF (bsize .lt. 1.E9_dp2) THEN
          ibsize = bsize/1.E6_dp2
          label = " Mb"
      ELSE
          ibsize = bsize/1.E9_dp2
          label = " Gb"
      END IF

      DO m = 1, 2
         IF (m .eq. 1) n = 6
         IF (m .eq. 2) n = nthreed
         WRITE (n, '(2x,a,i4,a,i3,a)') 'Block dim: ', mblk,
     1         '^2  Preconditioner size: ', ibsize, TRIM(label)
      END DO
!
!     COMPUTE AND STORE BLOCKS (MN X MN) FOR PRECONDITIONER
!
      CALL second0(time_on)

      IF (ALLOCATED(rzl_save)) DEALLOCATE (rzl_save)

      ALLOCATE (gc_save(ns,0:ntor,0:mpol1,ntyptot), stat=istat)
      IF (istat .ne. 0)
     1    STOP 'Allocation error: gc_save in compute_blocks'

      IF (ALLOCATED(block_diag))
     1    DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
      IF (istat .ne. 0) STOP 'Deallocation error in compute blocks'
      ALLOCATE (
     1 block_diag(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     2 block_plus(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     3 block_mins(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     4          stat=istat)
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'Allocation error(1) = ', istat,
     1               ': blocks in compute_blocks'
      END IF

      IF (ALLOCATED(lu_save))
     1   DEALLOCATE(lu_save, lv_save, r1_save, ru_save,
     2              rv_save, rcon_save, z1_save, zu_save,
     3              zv_save, zcon_save)
      ALLOCATE(
     1   lu_save(ns*nzeta,ntheta3,0:1), lv_save(ns*nzeta,ntheta3,0:1),
     1   ru_save(ns*nzeta,ntheta3,0:1), rv_save(ns*nzeta,ntheta3,0:1),
     2   r1_save(ns*nzeta,ntheta3,0:1), rcon_save(ns*nzeta,ntheta3,0:1),
     3   zu_save(ns*nzeta,ntheta3,0:1), zv_save(ns*nzeta,ntheta3,0:1),
     4   z1_save(ns*nzeta,ntheta3,0:1), zcon_save(ns*nzeta,ntheta3,0:1),
     5   stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error(2) in compute_blocks!'

      block_plus = 0; block_mins = 0
      IF (.not.lfreeb) block_diag(:,:,:,:,:,:,ns) = 0

      iequi_save = iequi;  iequi = 2                     !Signals funct3d NOT to precondition forces
      ncurr_save = ncurr
      IF (.not.liota) ncurr = 0                          !iota NOT pushed in time; should not contribute to Hessian

      nmax_jog = 2*ntmax
      IF (lasym) nmax_jog = ntyptot
      IF (lasym) lsweep_fast = .false.      !Not yet implemented for asym
      IF (ljog_test) THEN
         nmax_jog = ntyptot
         IF (lsweep_fast) THEN
            PRINT *,' Changing lsweep_fast to FALSE for jog test'
            lsweep_fast = .false.
         END IF
      END IF
!
!     GENERAL (SLOW) METHOD: ASSUMES NO SYMMETRIES OF R, Z COEFFICIENTS
!
      IF (.not.lsweep_fast)
     1   CALL sweep3_blocks (xc, xcdot, gc, gc_save, nmax_jog)

!
!     SYMMETRIC (R,Z, FAST) METHOD
!     Note: to compare sweep2 with sweep3 results,uncomment (1) "lsweep_fast=true" and
!           (2) "ljog_test =" lines in SWEEP2_BLOCKS, and (3) set block_ => blockl_ in SWEEP2_BLOCK
!      lsweep_fast = .true.
      IF (lsweep_fast)
     1   CALL sweep2_blocks (xc, xcdot, gc, gc_save, nmax_jog)

      DEALLOCATE(lu_save, lv_save, r1_save, ru_save, rv_save,
     1   rcon_save, z1_save, zu_save, zv_save, zcon_save)

      IF (lforbal) THEN
          block_diag(0,0,lamtype,:,:,:,:) = 0
          block_plus(0,0,lamtype,:,:,:,:) = 0
          block_mins(0,0,lamtype,:,:,:,:) = 0
      END IF

      CALL Compare_Jog_to_Analytic(xc)

!
!     FIXED BOUNDARY: R,Z NOT VARIED AT EDGE
!
      IF (.not.lfreeb) THEN
         block_mins(:,:,1:rztype,:,:,:,ns) = 0
         block_plus(:,:,:,:,:,1:rztype,ns-1) = 0
         block_diag(:,:,1:rztype,:,:,lamtype:,ns) = 0
      END IF

!
!     FILL IN NON-ZERO DIAGONAL ELEMENTS FOR (M=0,N), (M,N=0)
!     SIN MODES AND js=1, M > 0 MODES
!
      diag_val = block_diag(0,0,rcc,0,0,rcc,:)
      diag_lam = block_diag(0,1,zsc+2*ntmax,0,1,zsc+2*ntmax,:)
      diag_lam(1) = diag_lam(2)

!
!     M=1 constraint in RESIDUE: rsc(m=1) = zcc(m=1)
!     MAKE rsc CARRY PERTURBATION, AND SET zcc(m=1) = 0 BELOW
!     IT IS SET = RSC IN RESIDUE
!
      IF (lasym) THEN
         block_diag(:,:,:,:,1,rsc,:) = block_diag(:,:,:,:,1,rsc,:)
     1                               + block_diag(:,:,:,:,1,ntmax+zcc,:)
         block_mins(:,:,:,:,1,rsc,:) = block_mins(:,:,:,:,1,rsc,:)
     1                               + block_mins(:,:,:,:,1,ntmax+zcc,:)
         block_plus(:,:,:,:,1,rsc,:) = block_plus(:,:,:,:,1,rsc,:)
     1                               + block_plus(:,:,:,:,1,ntmax+zcc,:)
      END IF

      DO n = 0, ntor
         block_diag(n,0,zsc+ntmax,n,0,zsc+ntmax,:) = diag_val        !zsc(m=0,n)
         IF (.not.liota .or. n.ne.0)
     1   block_diag(n,0,zsc+2*ntmax,n,0,zsc+2*ntmax,:) = diag_lam    !lsc(m=0,n): n=0 stores iota if liota=.true.
         IF (.not.lasym) CYCLE
         block_diag(n,0,rsc,n,0,rsc,:) = diag_val
         block_diag(n,1,zcc+ntmax,n,1,zcc+ntmax,:) = diag_val        !zcc(m=1,n=0) = 0 (= rsc in residue)
         IF (jlam(0).gt.1)                    !lambda (m=0, js=1) not pushed, extrapolated: force = 0
     1   block_diag(n,0,zcc+2*ntmax,n,0,zcc+2*ntmax,1) = diag_lam(1)
      END DO

      IF (lasym) THEN
         block_diag(0,0,zcc+2*ntmax,0,0,zcc+2*ntmax,:) = diag_lam    !lcc(m=0,n=0)
      END IF

      IF (lthreed) THEN
         DO n = 0, ntor
            block_diag(n,0,rss,n,0,rss,:) = diag_val
            IF (jlam(0) .gt. 1)                         !lambda (m=0, js=1) not pushed...
     1      block_diag(n,0,zcs+2*ntmax,n,0,zcs+2*ntmax,1) =  diag_lam(1)
            IF (.not.lasym) CYCLE
            block_diag(n,0,zss+ntmax,n,0,zss+ntmax,:) = diag_val
            block_diag(n,0,zss+2*ntmax,n,0,zss+2*ntmax,:) =  diag_lam
         END DO
         DO m = 0, mpol1
            block_diag(0,m,rss,0,m,rss,:) = diag_val
            block_diag(0,m,zcs+ntmax,0,m,zcs+ntmax,:) = diag_val
            block_diag(0,m,zcs+2*ntmax,0,m,zcs+2*ntmax,:) =  diag_lam
            IF (.not.lasym) CYCLE
            block_diag(0,m,rcs,0,m,rcs,:) = diag_val
            block_diag(0,m,zss+ntmax,0,m,zss+ntmax,:) = diag_val
            block_diag(0,m,zss+2*ntmax,0,m,zss+2*ntmax,:) = diag_lam
         END DO
!
!        M=1 zcs force constraint: f(zcs) = 0, for all n, m=1
!        Also, zero variations with respect to zcs for all n, m=1
!
         block_diag(1:ntor,1,zcs+ntmax,:,:,:,:) = 0
         block_plus(1:ntor,1,zcs+ntmax,:,:,:,:) = 0
         block_mins(1:ntor,1,zcs+ntmax,:,:,:,:) = 0

         block_diag(:,:,:,1:ntor,1,zcs+ntmax,:) = 0
         block_plus(:,:,:,1:ntor,1,zcs+ntmax,:) = 0
         block_mins(:,:,:,1:ntor,1,zcs+ntmax,:) = 0

         DO n = 1, ntor
            block_diag(n,1,zcs+ntmax,n,1,zcs+ntmax,:) = diag_val
         END DO

      END IF

!
!     EDGE PEDESTAL FOR FREE-BOUNDARY
!
      IF (lfreeb) THEN
!         block_diag(:,:,1:rztype,:,:,1:rztype,ns) =
!     1   block_diag(:,:,1:rztype,:,:,1:rztype,ns)*1.05_dp2
         DO ntype = 1, rztype
            DO m = 0, mpol1
               DO n = 0, ntor
                  block_diag(n,m,ntype,n,m,ntype,ns) =
     1            block_diag(n,m,ntype,n,m,ntype,ns)*1.02_dp2
               END DO
            END DO
         END DO
      END IF

!
!     CHECK: M >= 1 MODES ARE NOT EVOLVED AT ORIGIN (js = 1)
!
      IF (ANY(block_plus(:,1:,:,:,:,:,1) .ne. zero)) THEN
          WRITE (6,*) 'block_plus(1) != 0'
          block_plus(:,1:,:,:,:,:,1) = 0
      END IF
      IF (ANY(block_mins(:,:,:,:,1:,:,2) .ne. zero)) THEN
          WRITE (6,*) 'block_mins(2) != 0'
          block_mins(:,:,:,:,1:,:,2) = 0
      END IF
!
!     PUT FINITE VALUES ON DIAGONAL AT JS=1, FOR M>=1, TO AVOID SINGULAR HESSIAN
!
      DO m = 1, mpol1
         DO n = 0, ntor
            DO ntype = 1, ntyptot
               block_diag(n,m,ntype,n,m,ntype,1) =
     1         block_diag(n,m,ntype,n,m,ntype,2)
            END DO
         END DO
      END DO

      IF (lthreed .and. jlam(0).eq.1) THEN
!        lcs m=0 mode pushed at origin;
!        must eliminate (near) zero eigenvalue due to only "half" mesh
!        value of lambda(n,m=0) really determined, so don't have an independent
!        force for origin value; pick small pedestal factor to avoid singular Hessian
         ntype = zcs+2*ntmax
         block_diag(:,0,ntype,:,0,ntype,1) = 1.01_dp2*
     1   block_diag(:,0,ntype,:,0,ntype,1)
      END IF

      IF (liota) THEN
         ntype = zsc + 2*ntmax
         block_diag(0,0,ntype,0,0,ntype,1) = diag_val(1)     !Not evolved in time...
      END IF

      CALL second0(time_off)
      DO m = 1, 2
         IF (m .eq. 1) n = 6
         IF (m .eq. 2) n = nthreed
         WRITE (n,'(1x,a,f10.2,a)')' Time to compute blocks: ',
     1      time_off - time_on,' s'
      END DO

!
!     CLEANUP
!
      iequi = iequi_save; ncurr = ncurr_save

!
!     Compute factorization of Hessian
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk, stat=ntype)
      ALLOCATE (ipiv_blk(mblk,ns), stat=ntype)
      IF (ntype .ne. 0) STOP 'Allocation error2 in block_precond'

      IF (l_backslv) THEN
!        Save blocks for checking in l_backslv=true loop below
         ALLOCATE (
     1   block_dsave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     2   block_msave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     3   block_psave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     4   stat = istat)
         IF (istat .ne. 0) THEN
            WRITE (6,*) 'Allocation error in l_backslv block: stat = ',
     1                istat
            l_backslv = .false.
         ELSE
            block_dsave = block_diag;  block_msave = block_mins
            block_psave = block_plus
         END IF
      END IF

      CALL second0(time_on)
      CALL blk3d_factor(block_diag, block_mins, block_plus,
     1                  ipiv_blk, mblk, ns)
      CALL second0(time_off)
      DO m = 1, 2
          IF (m .eq. 1) n = 6
          IF (m .eq. 2) n = nthreed
          WRITE (n,'(1x,a,f10.2,a)')' Time to factor blocks:  ',
     1                       time_off - time_on,' s'
      END DO

      IF (.not.l_backslv) DEALLOCATE (gc_save)

      CALL second0(tprec2doff)

      timer(tprec2d) = timer(tprec2d) + (tprec2doff - tprec2don)

      END SUBROUTINE compute_blocks


      SUBROUTINE sweep3_blocks (xc, xcdot, gc, gc_save, nmax_jog)
      USE vmec_main, ONLY: r01, z01
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: xc, xcdot,
     1                                                     gc, gc_save
      INTEGER :: nmax_jog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec2) :: eps, hj
      INTEGER :: radial_pts((ns+2)/3)
      INTEGER :: js, i, istat, nsize, mesh, lamtype, rztype, lfor
      LOGICAL, PARAMETER :: lscreen = .false.
C-----------------------------------------------
!
!     COMPUTE EVERY 3rd RADIAL POINT FOR EACH MESH STARTING AT js=1,2,3, RESPECTIVELY
!

      WRITE (6, *)
     1   " Using non-symmetric sweep to compute Hessian elements"

      eps = SQRT(EPSILON(eps))/3
      rztype = 2*ntmax
      lamtype = rztype+1
!     indicates forces to be computed using jogs...; others computed analytically
      lfor = nmax_jog

      xcdot = 0
      n_2d = 0;  m_2d = 0
      loglam_blks  = .not.lasym         !!Don't have lambda working yet for asymmetric case

!
!     CALL FUNCT3D FIRST TO STORE INITIAL UN-PRECONDITIONED FORCES
!     THIS WILL CALL LAMBLKS (loglam_blks = .true.) WHICH RESETS loglam_blks = .false.
!
      CALL funct3d (lscreen, istat)
      gc_save = gc

      MESH_LOOP: DO mesh = 1, 3

      radial_pts(1) = mesh
      nsize = 1

      DO js = 2, SIZE(radial_pts)
         radial_pts(js) = radial_pts(js-1) + 3
         IF (radial_pts(js) .gt. ns) EXIT
         nsize = nsize + 1
      END DO


!
!     PERFORM "JOGS" FOR ALL VARIABLES EVERY 3rd RADIAL POINT FROM EDGE
!     FOR Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss
!     POSSIBLY, ONLY NEED TO DO THIS FOR m <= 3 - 4, EXTRAPOLATE ABOVE THIS ?
!

!     Lambda jogs are implemented analytically
!     Except for lasym case

      istat = 0

      VAR2D_TYPE: DO ntype_2d = 1, nmax_jog
         IF (ntype_2d .lt. lamtype) THEN
            hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
         ELSE
            hj = eps
         END IF
         N2D: DO n_2d = 0, ntor
            M2D: DO m_2d = 0, mpol1
               DO i = 1, nsize
                  js = radial_pts(i)
                  xc(js,n_2d,m_2d,ntype_2d) =
     1            xc(js,n_2d,m_2d,ntype_2d) + hj
                  xcdot(js,n_2d,m_2d,ntype_2d) = hj
               END DO

               CALL funct3d (lscreen, istat)
               IF (istat .ne. 0) STOP 'Error computing Hessian jog!'

!
!              COMPUTE PRECONDITION (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF FORM (FIXED mn FOR SIMPLICITY):
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == block_mins; dj == block_diag; bj = block_plus
!
!              THUS:  d(js) = dF(js)/hj(js)
!                     b(js-1) = dF(js-1)/hj(js)
!                     a(js+1) = dF(js+1)/hj(js)
!
               DO i = 1, nsize
                  js = radial_pts(i)
                  xc(js,n_2d,m_2d,ntype_2d) =
     1            xc(js,n_2d,m_2d,ntype_2d) - hj
                  xcdot(js,n_2d,m_2d,ntype_2d) = 0

                  IF (js.eq.ns .and. .not.lfreeb) THEN
                     IF (ntype_2d .lt. lamtype) THEN
!                       PUT NON-ZERO DIAGONAL ELEMENT TO AVOID SINGULARITY
                        IF (ns .gt. 3)
     1            block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns) =
     2            block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns-3)
                        IF (ns .le. 3)
     1            block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns)
     2                   = 1
                     ELSE
                        block_diag(:,:,lamtype:,n_2d,m_2d,ntype_2d,ns) =
     1                  (gc(ns,:,:,lamtype:) -
     2                   gc_save(ns,:,:,lamtype:))/hj
                     END IF
                  ELSE
                     block_diag(:,:,1:lfor,n_2d,m_2d,ntype_2d,js) =
     1               (gc(js,:,:,1:lfor) - gc_save(js,:,:,1:lfor))/hj
                  END IF
!
!                 FOR OFF-DIAGONAL ELEMENTS, NEED TO ADJUST js INDICES +/- 1
!
                  IF (js .gt. 1)
     1               block_plus(:,:,1:lfor,n_2d,m_2d,ntype_2d,js-1) =
     2           (gc(js-1,:,:,1:lfor) - gc_save(js-1,:,:,1:lfor))/hj
                  IF (js .lt. ns)
     1               block_mins(:,:,1:lfor,n_2d,m_2d,ntype_2d,js+1) =
     2           (gc(js+1,:,:,1:lfor) - gc_save(js+1,:,:,1:lfor))/hj
               END DO
            END DO M2D
         END DO N2D
      END DO VAR2D_TYPE

      END DO MESH_LOOP

      END SUBROUTINE sweep3_blocks

      SUBROUTINE sweep2_blocks (xc, xcdot, gc, gc_save, nmax_jog)
      USE vmec_main, ONLY: r01, z01
      USE fbal
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: xc, xcdot,
     1                                                     gc, gc_save
      INTEGER :: nmax_jog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec2) :: eps, hj
      INTEGER :: radial_pts((ns+1)/2), js
      INTEGER :: i, istat, nsize, mesh, lamtype, rztype, lfor
      LOGICAL, PARAMETER :: lscreen = .false.
C-----------------------------------------------
!
!     COMPUTE EVERY 2nd RADIAL POINT FOR EACH MESH STARTING AT js=2,1 RESPECTIVELY
!     USES SYMMETRY OF HESSIAN WITH RESPECT TO R AND Z
!
!     TO CHECK THESE VALUE, SET BLOCK_ TO BLOCKL_ AND RUN BOTH THIS ROUTINE
!     AND THE SWEEP3_BLOCKS ROUTINE

      WRITE (6, *)
     1   " Using symmetric (fast) sweep to compute Hessian elements"

      eps = SQRT(EPSILON(eps))/3   !Use eps = 10*eps if central difference used
      rztype = 2*ntmax
      lamtype = rztype+1
!     indicates forces to be computed using jogs...; others computed analytically
!     note: can NOT compute lambda forces here, since they are not symmetric
      lfor = rztype

      xcdot = 0
      n_2d = 0;  m_2d = 0
      loglam_blks  = .not.lasym         !!Don't have lambda working yet for asymmetric case

!
!     CALL FUNCT3D FIRST TO STORE INITIAL UN-PRECONDITIONED FORCES
!     THIS WILL CALL LAMBLKS (loglam_blks = .true.) WHICH RESETS loglam_blks = .false.
!
!     To compare SWEEP2 with SWEEP3, uncomment the two ljog_test lines below (needed
!     to allocate blockl_ arrays and store lamblks calculations). Also,
!     change all block_ -> blockl_ below

!      ljog_test = .true.
      CALL funct3d (lscreen, istat)
!      ljog_test = .false.

      gc_save = gc

      MESH_LOOP: DO mesh = 0, 1             !even/odd spatial radial nodes

      radial_pts(1) = 2-mesh
      nsize = 1

      DO js = 2, SIZE(radial_pts)
         radial_pts(js) = radial_pts(js-1) + 2
         IF (radial_pts(js) .le. ns) THEN
            nsize = nsize + 1
         ELSE
            EXIT
         END IF
      END DO

!
!     PERFORM "JOGS" FOR ALL VARIABLES EVERY 2nd RADIAL POINT FROM EDGE
!     FOR Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss
!     PROBABLY ONLY NEED TO DO THIS FOR m <= 3 - 4, EXTRAPOLATE ABOVE THIS...
!
!     NOTE: FOR THESE FORMULAE TO WORK, HJ MUST BE A CONSTANT
!           OTHERWISE, THE FORMULAE MUST BE AMENDED FOR VARIABLE HJ
!

!     Lambda jogs are implemented analytically
!     Except for lasym case

      istat = 0

      VAR2D_TYPE: DO ntype_2d = 1, nmax_jog
         IF (ntype_2d .lt. lamtype) THEN
            hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
         ELSE
            hj = eps
         END IF
         N2D: DO n_2d = 0, ntor
            M2D: DO m_2d = 0, mpol1
               DO i = 1, nsize
                  js = radial_pts(i)
                  xc(js,n_2d,m_2d,ntype_2d) =
     1            xc(js,n_2d,m_2d,ntype_2d) + hj
                  xcdot(js,n_2d,m_2d,ntype_2d) = hj
               END DO

               CALL funct3d (lscreen, istat)
               IF (istat .ne. 0) STOP 'Error computing Hessian jog!'

!
!              COMPUTE PRECONDITION (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF FORM (FIXED mn FOR SIMPLICITY):
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == block_mins; dj == block_diag; bj = block_plus
!
!              FOR EVEN mesh = 0 (js = 2, 4, 6....), perturbation = (0, 1, 0, 1, ....)
!
!              F(1) = b1
!              F(2) = d2
!              F(3) = a3 + b3
!              F(2k) = d(2k);   F(2k+1) = a(2k+1) + b(2k+1)
!
!              FOR ODD  mesh = 1 (js = 1, 3, 5, ...), perturbation = (1, 0, 1, 0, ...)
!
!              F(1) = d1
!              F(2) = a2 + b2
!              F(3) = d3
!              F(2k) = a(2k) + b(2k);  F(2k+1) = d(2k+1)
!
!              Use symmetry, b(k+1) = TRANSPOSE(a(k)) TO SOLVE THIS SYSTEM AS FOLLOWS:
!
!              d(k) = F(k)  for all k (this is easy)
!
!
!
               DO i = 1, nsize
                  js = radial_pts(i)
                  xc(js,n_2d,m_2d,ntype_2d) =
     1            xc(js,n_2d,m_2d,ntype_2d) - hj
                  xcdot(js,n_2d,m_2d,ntype_2d) = 0

                  IF (js.eq.ns .and. .not.lfreeb) THEN
                     IF (ntype_2d .lt. lamtype) THEN
!                       PUT NON-ZERO DIAGONAL ELEMENT TO AVOID SINGULARITY
                        IF (ns .gt. 3)
     1          block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns) =
     2          block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns-2)
                        IF (ns .le. 3)
     1            block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,ns)
     2                   = 1
                     ELSE
                      block_diag(:,:,lamtype:,n_2d,m_2d,ntype_2d,ns) =
     1                  (gc(ns,:,:,lamtype:) -
     2                   gc_save(ns,:,:,lamtype:))/hj
                     END IF
                  ELSE
                     block_diag(:,:,1:lfor,n_2d,m_2d,ntype_2d,js) =
     1               (gc(js,:,:,1:lfor) - gc_save(js,:,:,1:lfor))/hj
                  END IF
!
!                 STORE TEMPORARY OFF-DIAGONAL ELEMENTS IN block_plus
!                 (actual a,b will be computed and stored in symm_blocks)
!                 NOTE: boundary conditions:
!                     a(1) = 0; => b(1) = F(1);  [block_plus(1) == b(1)]
!                     b(ns) = 0
!
!                 EVEN MESH: js = 2, 4, 6, ...
!                 block_plus(js-1) == a(js-1) + b(js-1) = F(js-1)
!
!                 ODD MESH:  js = 1, 3, 5, ...
!                 block_plus(js+1) == a(js+1) + b(js+1) = F(js+1)
!
                  IF (js.gt.1 .and. mesh.eq.0) THEN               !Even mesh, starts at js=2
                     block_plus(:,:,1:lfor,n_2d,m_2d,ntype_2d,js-1) =
     1           (gc(js-1,:,:,1:lfor) - gc_save(js-1,:,:,1:lfor))/hj
                  ELSE IF (js.lt.ns .and. mesh.eq.1) THEN         !Odd mesh, starts at js=1
                     block_plus(:,:,1:lfor,n_2d,m_2d,ntype_2d,js+1) =
     2           (gc(js+1,:,:,1:lfor) - gc_save(js+1,:,:,1:lfor))/hj
                  END IF

               END DO
            END DO M2D
         END DO N2D
      END DO VAR2D_TYPE

      END DO MESH_LOOP

!
!     COMPUTE OFF-DIAGONAL BLOCKS USING FORMULAE GIVEN ABOVE
!     FOR js > 1:
!
!     b(js) = F(js) - a(js) == F(js) - TRANSPOSE(b(js-1))
!
      CALL symm_blocks (block_plus, block_mins, mnsize, rztype)

!
!     block_ (at lamtype) used to store Hessian of EQUIF
!
      IF (lforbal) THEN
         DO js = 2, ns1
         block_diag(0,1,rcc,:,:,1:rztype,js) =
     1      frcc_fac(js)*block_diag(0,1,rcc,:,:,1:rztype,js)
     2    + fzsc_fac(js)*block_diag(0,1,zsc+ntmax,:,:,1:rztype,js)

         block_diag(0,1,zsc+ntmax,:,:,1:rztype,js) =
     1      rru_fac(js)*(block_diag(0,0,lamtype,:,:,1:rztype,js)/4
     2                 - block_diag(0,1,rcc,:,:,1:rztype,js))

         block_diag(0,1,rcc,:,:,1:rztype,js) =
     1      rzu_fac(js)*(block_diag(0,0,lamtype,:,:,1:rztype,js)/4
     2                 + block_diag(0,1,rcc,:,:,1:rztype,js))

         block_plus(0,1,rcc,:,:,1:rztype,js) =
     1      frcc_fac(js)*block_plus(0,1,rcc,:,:,1:rztype,js)
     2    + fzsc_fac(js)*block_plus(0,1,zsc+ntmax,:,:,1:rztype,js)

         block_plus(0,1,zsc+ntmax,:,:,1:rztype,js) =
     1      rru_fac(js)*(block_plus(0,0,lamtype,:,:,1:rztype,js)/4
     2                 - block_plus(0,1,rcc,:,:,1:rztype,js))

         block_plus(0,1,rcc,:,:,1:rztype,js) =
     1      rzu_fac(js)*(block_plus(0,0,lamtype,:,:,1:rztype,js)/4
     2                 + block_plus(0,1,rcc,:,:,1:rztype,js))

         block_mins(0,1,rcc,:,:,1:rztype,js) =
     1      frcc_fac(js)*block_mins(0,1,rcc,:,:,1:rztype,js)
     2    + fzsc_fac(js)*block_mins(0,1,zsc+ntmax,:,:,1:rztype,js)

         block_mins(0,1,zsc+ntmax,:,:,1:rztype,js) =
     1      rru_fac(js)*(block_mins(0,0,lamtype,:,:,1:rztype,js)/4
     2                 - block_mins(0,1,rcc,:,:,1:rztype,js))

         block_mins(0,1,rcc,:,:,1:rztype,js) =
     1      rzu_fac(js)*(block_mins(0,0,lamtype,:,:,1:rztype,js)/4
     2                 + block_mins(0,1,rcc,:,:,1:rztype,js))
         END DO
      END IF

!
!     IMPOSE M=1 CONSTRAINTS (FROM TOTZSP: R,Z(m=1,j=1) = R,Z(m=1,j=2))
!
      block_diag(:,:,1:lfor,:,1,1:lfor,2) =
     1     block_diag(:,:,1:lfor,:,1,1:lfor,2)
     2   + block_mins(:,:,1:lfor,:,1,1:lfor,2)

      block_plus(:,0,1:lfor,:,1,1:lfor,1) =
     1     block_plus(:,0,1:lfor,:,1,1:lfor,1)
     2   + block_diag(:,0,1:lfor,:,1,1:lfor,1)

      block_diag(:,:,1:lfor,:,1,1:lfor,1) = 0
      block_mins(:,:,1:lfor,:,1,1:lfor,2) = 0
      block_diag(:,1,1:lfor,:,:,1:lfor,1) = 0
      block_plus(:,1,1:lfor,:,:,1:lfor,1) = 0

      END SUBROUTINE sweep2_blocks

      SUBROUTINE symm_blocks (block_p, block_m, mblk, rztype)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, rztype
      REAL(sp), DIMENSION(mblk,ntyptot,mblk,ntyptot,ns) ::
     1                    block_p, block_m
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, nt, ntp, nsmax
C-----------------------------------------------
!
!     ON INPUT, block_p CONTAIN FORCE PERTURBATIONS FOR ALL js (even/odd meshes combined)
!     ON OUTPUT, block_p, block_m contain the (symmetric) plus/minus blocks
!
!     NOTE: js MUST be the outer loop (all nt,ntp type combos must be processed at a given node FIRST)
!
      nsmax = ns
      IF (.not.lfreeb) nsmax = ns-1
      DO js = 2, nsmax
         DO nt = 1, rztype
            DO ntp = 1, rztype
               block_m(:,nt,:,ntp,js) =
     1         TRANSPOSE(block_p(:,ntp,:,nt,js-1))
               IF (js .eq. ns) CYCLE
               block_p(:,nt,:,ntp,js) = block_p(:,nt,:,ntp,js)
     1                                - block_m(:,nt,:,ntp,js)
            END DO
         END DO
      END DO

      END SUBROUTINE symm_blocks

      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE laprec
!     USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec2 = 8
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
      REAL(sp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(inout) ::
     1                       a, bm1, bp1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
      REAL(sp), POINTER :: amat(:,:), bmat(:,:)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c-----------------------------------------------------------------------
c
c  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  a                   : diagonal block
c  bp1, bm1            : lower, upper blocks (see equation below)
c
c  OUTPUT
c  ipiv                : pivot elements for kth block
c  a                   : a-1 LU factor blocks
c  bm1                 : q = a-1 * bm1 matrix
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
c
c     1. Start from row N and solve for x(N) in terms of x(N-1):
c
c        x(N) = -q(N)*x(N-1) + r(N)
c
c        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
c
c        where a(N)[-1] is the inverse of a(N)
c
c     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
c
c        x(l) = -q(l)*x(l-1) + r(l), in general, where:
c
c        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
c
c        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
c
c        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
c
c     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
c
c        x(1) = r(1)
c
c     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK, la_... ROUTINES)
c
c     1. CALL la_getrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
c     2. CALL la_getrs:   With multiple (mblk) right-hand sides, to do block inversion
c                         operation, A X = B  (stores result in B; here B is a matrix)
c
!      ndisk = mblk

c  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) STOP 'Error opening scratch file in blk3d'

c  main loop. load and process (backwards) block-rows nblocks to 1.


      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(nblocks) = ablk(nblocks)[-1] * bml
!
         amat => a(:,:,k);  ipivot => ipiv(:,k)
         CALL la_getrf (mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
         IF (k .eq. 1) EXIT

         bmat => bm1(:,:,k)
         CALL la_getrs('n', mblk, mblk, amat, mblk, ipivot,
     1                 bmat, mblk, ier)

         IF (ier .ne. 0) GOTO 305

!         CALL wrdisk(iunit, ql, ndisk, incnow, ibuph, incbu, ier)
!         IF (ier .ne. 0) GOTO 302

!
!      Update effective diagonal "a" matrix
!
         k1 = k-1
         a(:,:,k1)  = a(:,:,k1) - MATMUL(bp1(:,:,k1), bm1(:,:,k))

      END DO BLOCKS

!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING SUM(BP1 * SOURCE)
!     IS MUCH FASTER (UP TO 3X FASTER ON LINUX) THAN MATMUL OPERATION
!
      DO k = 1, nblocks
         IF (k .ne. nblocks) bp1(:,:,k) = TRANSPOSE(bp1(:,:,k))
         IF (k .ne. 1)       bm1(:,:,k) = TRANSPOSE(bm1(:,:,k))
      END DO

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '
     1                        , k
      IF (ier < 0)WRITE (6,'(i4, a)')
     1    ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)')
     1    ier, 'th diagonal factor exactly zero'
      STOP
  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

c  destroy disk file and return. ---------------------------------------

  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_factor

      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source,
     1                     ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      USE laprec
!      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec2 = 8
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(sp), DIMENSION(mblk,mblk,nblocks), INTENT(in) ::
     1                       bp1, qblk
      REAL(sp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(in) ::
     1                       ablk
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(inout)
     1                  :: source
      INTEGER, POINTER  :: ipivot(:)
      REAL(sp), POINTER :: amat(:,:)
      REAL(sp) :: source_sp(mblk)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c-----------------------------------------------------------------------
c
c  this subroutine solves a block-tridiagonal system of equations, using
c  the ABLK, QBLK factors from blk3d_factor,
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  bp1                 : upper blocks (see equation below)
c  ipiv                : pivot elements for kth block
c  ablk                : a-1 blocks
c  qblk                : q = a-1 * bm1
c  source              : input right side
c
c  OUTPUT
c  source              : Solution x of A x = source
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
c
c     1. Start from row N and solve for x(N) in terms of x(N-1):
c
c        x(N) = -q(N)*x(N-1) + r(N)
c
c        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
c
c        where a(N)[-1] is the inverse of a(N)
c
c     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
c
c        x(l) = -q(l)*x(l-1) + r(l), in general, where:
c
c        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
c
c        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
c
c        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
c
c     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
c
c        x(1) = r(1)
c
c     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK, la_... ROUTINES)
c
c     1. CALL la_getrs:   With single right hand side (source) to solve A x = b (b a vector)
c                         Faster than la_gesl
!      ndisk = mblk

c  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) STOP 'Error opening scratch file in blk3d'

c  main loop. load and process (backwards) block-rows nblocks to 1.
!  note: about equal time is spent in calling la_getrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_sp(:) = source(:,k)
         ipivot => ipiv(:,k);   amat => ablk(:,:,k)
         CALL la_getrs('n', mblk, 1, amat, mblk,
     1                    ipivot, source_sp, mblk, ier)
         source(:,k) = source_sp(:)

         IF (ier .ne. 0) GOTO 305
         IF (k .eq. 1) EXIT

!
!        SUM IS MUCH FASTER THAN MATMUL OPERATION; NOTE WE TRANSPOSED
!        BP1 (AND BM1) IN BLK3D_FACTOR TO MAKE FIRST INDEX FASTEST VARYING
!        source(:,k-1) = source(:,k-1) - MATMUL(bp1(:,:,k-1),source(:,k))
!
         DO k1 = 1,mblk
            source(k1,k-1) = source(k1,k-1)
     1                     - SUM(bp1(:,k1,k-1)*source(:,k))
         END DO
      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks
!         CALL rddisk (iunit, ql, ndisk, incnow, ibuph, ier)
!         IF (ier .ne. 0) GOTO 303
!         ibuph = ibuph - incbu
!
!        NOTE: qblk in loop is transpose of qblk in matmul; faster this way
!        source(:,k) = source(:,k) - MATMUL(qblk(:,:,k),source(:,k-1))
!
         DO k1 = 1,mblk
            source(k1,k) = source(k1,k)
     1                   - SUM(qblk(:,k1,k)*source(:,k-1))
         END DO

      END DO

      GOTO 400

c  error returns. ------------------------------------------------------

  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

c  destroy disk file and return. ---------------------------------------

  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_slv

      SUBROUTINE Compare_Jog_to_Analytic (xc)
      IMPLICIT NONE
      REAL(rprec2),DIMENSION(ns,0:ntor,0:mpol1,3*ntmax) :: xc
      REAL(sp), PARAMETER :: tol1 = 1.E-2, tol2 = 1.E-3
      INTEGER :: ntype, istat, n, m, np, mp, js
      REAL(sp) :: t1, t2, t3
      LOGICAL  :: lcond
!
!     THIS ROUTINE PRINTS OUT A COMARISON OF THE 'JOGGED' HESSIAN
!     (BLOCK_...) WITH THE ANALYTICALLY COMPUTED VALUES, BLOCKL_...
!
      IF (.not.ALLOCATED(blockl_diag)) RETURN

      WRITE (6, *) ' Hessian comparison stored in FORT.35 file'
      WRITE (35, '(3x,a,/)')'LEGEND: FORCE (UNPRIMED)  VARIABLE (PRIME)'
      WRITE (35, *)'   TOTAL NUMBER OF R,Z HESSIAN ELEMENTS: ',
     1          (2*ntmax*mnsize)**2 * (3*ns-2)
      WRITE (35, '(/,2a)')'   js   N   M   N''  M'' TYP TYP''',
     1               '  BLOCK(num)    BLOCK(ana)    XC(N'',M'',TYPE'')'

      DO ntype = 1, 3*ntmax
      DO istat = 1, 3*ntmax
!      DO ntype = 1, 2*ntmax
!      DO istat = 1, 2*ntmax

!         IF (istat .ne. ntype) CYCLE    !ONLY DIAG ELEMENTS FOR NOW
!         GOTO 2000

         DO n = 0, ntor
            DO m = 0, mpol1
!               IF (m.ne.1 .or. n.ne.0) CYCLE   !CHECK LFORBAL ONLY
               DO np = 0, ntor
                  DO mp = 0,mpol1
                     lcond = (ntype.gt.2*ntmax .or. istat.gt.2*ntmax)  !Lambda force/pert
                     IF (lcond) CYCLE
!                     IF (.not.lcond) CYCLE
                     DO js = 1, ns-1
                        t1 = block_mins(n,m,ntype,np,mp,istat,js)
                        t2 = blockl_mins(n,m,ntype,np,mp,istat,js)
                        t3 = xc(js,np,mp,istat)
                        IF ((ABS(t1 - t2)/(ABS(t1)+ABS(t2)+1.E-8)
     1                     .gt. tol1) .and. (ABS(t1)+ABS(t2)).gt.tol2
     1                     .and. js.gt.1)
     1                  WRITE (35, 1234) '-', js, n, m, np, mp,
     1                                   ntype, istat, t1, t2, t3
                        t1 = block_diag(n,m,ntype,np,mp,istat,js)
                        t2 = blockl_diag(n,m,ntype,np,mp,istat,js)
                        IF ((ABS(t1 - t2)/(ABS(t1)+ABS(t2)+1.E-8)
     1                     .gt. tol1) .and. (ABS(t1)+ABS(t2)).gt.tol2)
     1                  WRITE (35, 1234) 'd',js, n, m, np, mp, ntype,
     1                                   istat, t1, t2, t3
                        t1 = block_plus(n,m,ntype,np,mp,istat,js)
                        t2 = blockl_plus(n,m,ntype,np,mp,istat,js)
                        IF ((ABS(t1 - t2)/(ABS(t1)+ABS(t2)+1.E-8)
     1                  .gt. tol1) .and. (ABS(t1)+ABS(t2)).gt.tol2
     1                  .and. js.lt.ns)
     1                  WRITE (35, 1234) '+', js, n, m, np, mp, ntype,
     1                                   istat, t1, t2, t3
                     END DO
                  END DO
               END DO
            END DO
         END DO

         CYCLE

! 2000    block_mins(:,:,ntype,:,:,istat,1:ns) =
!     1   blockl_mins(:,:,ntype,:,:,istat,1:ns)
!         block_plus(:,:,ntype,:,:,istat,1:ns) =
!     1   blockl_plus(:,:,ntype,:,:,istat,1:ns)
!         block_diag(:,:,ntype,:,:,istat,1:ns) =
!     1   blockl_diag(:,:,ntype,:,:,istat,1:ns)

      END DO
      END DO
 1234 FORMAT(a,7i4,1p,3e14.4)

      DEALLOCATE (blockl_diag, blockl_mins, blockl_plus)

      END SUBROUTINE Compare_Jog_to_Analytic

      END MODULE precon2d
