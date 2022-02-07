      SUBROUTINE wrout_cdf(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, 
     1   bsupu, rzl_array, gc_array, ier_flag, lwrite)
      USE vmec_main, p5 => cp5, two => c2p0
      USE vmec_params
      USE vmercier
      USE vmec_persistent
      USE vsvd
      USE vspline
      USE xstuff
      USE vmec_io
      USE realspace
      USE ezcdf
      USE read_wout_mod, ONLY: vn_version, vn_extension, vn_mgrid,
     1  vn_magen, vn_therm, vn_gam, vn_maxr, vn_minr, vn_maxz, vn_fp,
     2  vn_radnod, vn_polmod, vn_tormod, vn_maxmod, vn_maxit, vn_actit,
     3  vn_asym, vn_recon, vn_free, vn_error, vn_aspect, vn_beta, 
     4  vn_pbeta, vn_tbeta, vn_abeta, vn_b0, vn_rbt0, vn_maxmod_nyq,
     5  vn_rbt1, vn_sgs, vn_lar, vn_modB, vn_ctor, vn_amin, vn_Rmaj, 
     6  vn_vol, vn_mse, vn_thom, 
     7  vn_pmod, vn_tmod, vn_pmod_nyq, vn_tmod_nyq,
     7  vn_racc, vn_zacs, vn_racs, vn_zacc, vn_iotaf,
     8  vn_presf, vn_phi, vn_phipf, vn_jcuru, vn_jcurv, vn_iotah,
     9  vn_mass, vn_presh, vn_betah, vn_buco, vn_bvco, vn_vp, vn_specw, 
     A  vn_phip, vn_jdotb, vn_overr, vn_bgrv, vn_merc, vn_mshear,
     B  vn_mwell, vn_mcurr, vn_mgeo, vn_equif, vn_fsq, vn_wdot, 
     C  vn_extcur, vn_curlab, vn_rmnc, vn_zmns, vn_lmns, vn_gmnc, 
     D  vn_bmnc, vn_bsubumnc, vn_bsubvmnc, vn_bsubsmns, 
     E  vn_bsupumnc, vn_bsupvmnc, vn_rmns, vn_zmnc, vn_lmnc, vn_gmns,
     F  vn_bmns, vn_bsubumns, vn_bsubvmns, vn_bsubsmnc, vn_bsupumns, 
     G  vn_bsupvmns, vn_rbc, vn_zbs, vn_rbs, vn_zbc,
     H  ln_version, ln_extension, ln_mgrid,
     1  ln_magen, ln_therm, ln_gam, ln_maxr, ln_minr, ln_maxz, ln_fp,
     2  ln_radnod, ln_polmod, ln_tormod, ln_maxmod, ln_maxit, ln_actit,
     3  ln_asym, ln_recon, ln_free, ln_error, ln_aspect, ln_beta, 
     4  ln_pbeta, ln_tbeta, ln_abeta, ln_b0, ln_rbt0, ln_maxmod_nyq,
     5  ln_rbt1, ln_sgs, ln_lar, ln_modB, ln_ctor, ln_amin, ln_Rmaj, 
     6  ln_mse, ln_thom, ln_flp, ln_nobd, ln_nbset, ln_next, ln_nbfld,
     7  ln_pmod, ln_tmod, ln_pmod_nyq, ln_tmod_nyq, ln_racc, ln_zacs, 
     7  ln_racs, ln_zacc, ln_iotaf,
     8  ln_presf, ln_phi, ln_phipf, ln_jcuru, ln_jcurv, ln_iotah,
     9  ln_mass, ln_presh, ln_betah, ln_buco, ln_bvco, ln_vp, ln_specw, 
     A  ln_vol, ln_phip, ln_jdotb, ln_bgrv, ln_merc, ln_mshear,
     B  ln_mwell, ln_mcurr, ln_mgeo, ln_equif, ln_fsq, ln_wdot, 
     C  ln_extcur, ln_curlab, ln_rmnc, ln_zmns, ln_lmns, ln_gmnc, 
     D  ln_bmnc, ln_bsubumnc, ln_bsubvmnc, ln_bsubsmns, 
     E  ln_bsupumnc, ln_bsupvmnc, ln_rmns, ln_zmnc, ln_lmnc, ln_gmns,
     F  ln_bmns, ln_bsubumns, ln_bsubvmns, ln_bsubsmnc, ln_bsupumns, 
     G  ln_bsupvmns, ln_rbc, ln_zbs, ln_rbs, ln_zbc
      USE mgrid_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
      REAL(rprec), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)),           !reverse ns, mnmax for backwards compatibility
     1   INTENT(inout), TARGET :: rzl_array, gc_array
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu
      LOGICAL :: lwrite
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      CHARACTER*(*), PARAMETER, DIMENSION(1) ::
     1             r1dim = (/'radius'/), mn1dim = (/'mn_mode'/),
     2             mn2dim = (/'mn_mode_nyq'/),
     3             currg = (/'current_group'/)
      CHARACTER*(*), DIMENSION(2), PARAMETER :: 
     1             r2dim = (/'mn_mode','radius '/),
     1             r3dim = (/'mn_mode_nyq','radius     '/),
     2             currl = (/'clabel','curgrp'/)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, js, mn, lk, lr, 
     1           m, n, k, iwout0, n1, nwout
      REAL(rprec) :: dmult, tcosi, tsini, vversion, sgn, tmult
      REAL(rprec), POINTER, DIMENSION(:,:) :: rmnc, rmns, zmns, 
     1   zmnc, lmns, lmnc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: 
     1   gmnc, bmnc, gmns, bmns, 
     2   bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc
      REAL(rprec), DIMENSION(nznt) :: bmod
      REAL(rprec), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,
     1   rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
      REAL(rprec), DIMENSION(mnmax_nyq) :: gmn, bmn,
     1   bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
      CHARACTER*120 :: wout_file
!     ELIMINATE THESE EVENTUALLY
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: 
     1   bsupumnc, bsupumns, bsupvmnc, bsupvmns
C-----------------------------------------------
!
!  Pointer assignments for storage arrays
!
      n1 = MAX(1,ntmax/2)
      rmnc => rzl_array(:,:,1)            !!store COS(mu-nv) components
      zmns => rzl_array(:,:,1+n1)         !!store SIN(mu-nv)
      lmns => rzl_array(:,:,1+2*n1)       !!store SIN(mu-nv)

      IF (lasym) THEN
         rmns => gc_array(:,:,1)            !!store SIN(mu-nv)
         zmnc => gc_array(:,:,1+n1)         !!store COS(mu-nv)
         lmnc => gc_array(:,:,1+2*n1)       !!store COS(mu-nv)
      END IF

!
!     THIS SUBROUTINE CREATES THE NETCDF FILE WOUT WHICH
!     CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM
!     INTERNAL full REPRESENTATION)
!
!     IZETA (FULL), BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF)
!

      wout_file = 'wout_' // TRIM(input_extension) // '.nc'

      call cdf_open(nwout,wout_file,'w',iwout0)
      if (iwout0 .ne. 0) STOP 'Error opening WOUT.nc file VMEC WROUT'

      IF (.not. lrecon) THEN
         itse = 0
         imse2 = 0
      END IF

!================================
! Define Variables
!================================
!  Scalars 
      wout_file = version_
      READ (wout_file,*) vversion
      CALL cdf_define(nwout, vn_version, vversion)
      CALL cdf_define(nwout, vn_extension, TRIM(input_extension))
      CALL cdf_define(nwout, vn_mgrid, TRIM(mgrid_file))
      CALL cdf_define(nwout, vn_magen, wb)
      CALL cdf_define(nwout, vn_therm, wp)
      CALL cdf_define(nwout, vn_gam, gamma)
      CALL cdf_define(nwout, vn_maxr, rmax_surf)
      CALL cdf_define(nwout, vn_minr, rmin_surf)
      CALL cdf_define(nwout, vn_maxz, zmax_surf)
      CALL cdf_define(nwout, vn_fp, nfp)
      CALL cdf_define(nwout, vn_radnod, ns)
      CALL cdf_define(nwout, vn_polmod, mpol)
      CALL cdf_define(nwout, vn_tormod, ntor)
      CALL cdf_define(nwout, vn_maxmod, mnmax)
      CALL cdf_define(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_define(nwout, vn_maxit, niter)
      CALL cdf_define(nwout, vn_actit, itfsq)
      CALL cdf_define(nwout, vn_asym, lasym)
      CALL cdf_define(nwout, vn_recon, lrecon)
      CALL cdf_define(nwout, vn_free, lfreeb)
      CALL cdf_define(nwout, vn_error, ier_flag)
      CALL cdf_define(nwout, vn_aspect, aspect)
      CALL cdf_define(nwout, vn_beta, betatot)
      CALL cdf_define(nwout, vn_pbeta, betapol)
      CALL cdf_define(nwout, vn_tbeta, betator)
      CALL cdf_define(nwout, vn_abeta, betaxis)
      CALL cdf_define(nwout, vn_b0, b0)
      CALL cdf_define(nwout, vn_rbt0, rbtor0)
      CALL cdf_define(nwout, vn_rbt1, rbtor)
      CALL cdf_define(nwout, vn_sgs, NINT(signgs))
      CALL cdf_define(nwout, vn_lar, IonLarmor)
      CALL cdf_define(nwout, vn_modB, volAvgB)
      CALL cdf_define(nwout, vn_ctor, ctor)
      CALL cdf_define(nwout, vn_amin, Aminor_p)
      CALL cdf_define(nwout, vn_Rmaj, Rmajor_p)
      CALL cdf_define(nwout, vn_vol, volume_p)

      IF (lrecon) THEN
         CALL cdf_define(nwout, vn_mse, imse2)
         CALL cdf_define(nwout, vn_thom, itse)
      END IF

      IF (nextcur .eq. 0) nextcur = COUNT(extcur .ne. zero)
      CALL cdf_define(nwout, vn_nextcur, nextcur)
      CALL cdf_define(nwout, vn_extcur,
     1        extcur(1:nextcur), dimname=currg)
      CALL cdf_define(nwout, vn_mgmode, mgrid_mode)
      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_flp, nobser)
         CALL cdf_define(nwout, vn_nobd, nobd)
         CALL cdf_define(nwout, vn_nbset, nbsets)
         IF (nbsets .gt. 0) 
     1      CALL cdf_define(nwout,vn_nbfld,nbfld(1:nbsets))
      END IF

      IF (.not.lwrite) GO TO 800

! 1D Arrays

      CALL cdf_define(nwout, vn_pmod, xm, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_pmod, ln_pmod)
      CALL cdf_define(nwout, vn_tmod, xn, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_tmod, ln_tmod)
      CALL cdf_define(nwout, vn_pmod_nyq, xm_nyq, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_pmod_nyq, ln_pmod_nyq)
      CALL cdf_define(nwout, vn_tmod_nyq, xn_nyq, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_tmod_nyq, ln_tmod_nyq)

      CALL cdf_define(nwout, vn_racc, raxis_cc(0:ntor), 
     1                dimname=(/'n-tor'/))
      CALL cdf_setatt(nwout, vn_racc, ln_racc)
      CALL cdf_define(nwout, vn_zacs, zaxis_cs(0:ntor), 
     1                dimname=(/'n-tor'/))
      CALL cdf_setatt(nwout, vn_zacs, ln_zacs)
      IF (lasym) THEN
         CALL cdf_define(nwout, vn_racs, raxis_cs(0:ntor), 
     1                dimname=(/'n-tor'/))
         CALL cdf_setatt(nwout, vn_racs, ln_racs)
         CALL cdf_define(nwout, vn_zacc, zaxis_cc(0:ntor), 
     1                dimname=(/'n-tor'/))
         CALL cdf_setatt(nwout, vn_zacc, ln_zacc)
      END IF

      CALL cdf_define(nwout, vn_iotaf, iotaf(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_iotaf)
      CALL cdf_define(nwout, vn_presf, presf, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presf, ln_presf, units='Pa')
      CALL cdf_define(nwout, vn_phi, phi, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phi, ln_phi, units='wb')
      CALL cdf_define(nwout, vn_phipf, 
     1                phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phipf, ln_phipf)
      CALL cdf_define(nwout, vn_jcuru, 
     1                jcuru, dimname=r1dim)
      CALL cdf_define(nwout, vn_jcurv, 
     1                jcurv, dimname=r1dim)
 
      CALL cdf_define(nwout, vn_iotah, iotas(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotah, ln_iotah)
      CALL cdf_define(nwout, vn_mass, mass, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_mass, ln_mass)
      CALL cdf_define(nwout, vn_presh, pres(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presh, ln_presh, units='Pa')
      CALL cdf_define(nwout, vn_betah, beta_vol, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_buco, buco, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bvco, bvco, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_vp, vp(1:ns), 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_specw, specw, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_phip, 
     1                phips(1:ns), dimname=r1dim)
      CALL cdf_define(nwout, vn_overr, 
     2                overr(1:ns), dimname=r1dim)

      CALL cdf_define(nwout, vn_jdotb, jdotb,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bgrv, bdotgradv,
     1                dimname=r1dim)

      CALL cdf_define(nwout, vn_merc, Dmerc,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mshear, Dshear,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mwell, Dwell,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mcurr, Dcurr,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mgeo,
     1                Dgeod, dimname=r1dim)
      CALL cdf_define(nwout, vn_equif,
     1                equif, dimname=r1dim)

      CALL cdf_define(nwout, vn_fsq, fsqt(1:nstore_seq),
     1                dimname=(/'time'/))
      CALL cdf_define(nwout, vn_wdot, wdot(1:nstore_seq),
     1                dimname=(/'time'/))

      IF (lfreeb .and. nextcur.gt.0) THEN
         CALL cdf_define(nwout, vn_curlab,
     1        curlabel(1:nextcur), dimname=currl)
      ENDIF

! 2D Arrays
      ALLOCATE (gmnc(mnmax_nyq,ns), bmnc(mnmax_nyq,ns),
     1          bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns),
     2          bsubsmns(mnmax_nyq,ns), bsupumnc(mnmax_nyq,ns),
     3          bsupvmnc(mnmax_nyq,ns))
      CALL cdf_define(nwout, vn_rmnc, rmnc, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_rmnc, ln_rmnc, units='m')
      CALL cdf_define(nwout, vn_zmns, zmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_zmns, ln_zmns, units='m')
      CALL cdf_define(nwout, vn_lmns, lmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_lmns, ln_lmns)
      CALL cdf_define(nwout, vn_gmnc, gmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_gmnc, ln_gmnc)
      CALL cdf_define(nwout, vn_bmnc, bmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bmnc, ln_bmnc)
      CALL cdf_define(nwout, vn_bsubumnc, bsubumnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubumnc, ln_bsubumnc)
      CALL cdf_define(nwout, vn_bsubvmnc, bsubvmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubvmnc, ln_bsubvmnc)
      CALL cdf_define(nwout, vn_bsubsmns, bsubsmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubsmns, ln_bsubsmns)

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM - CAN COMPUTE FROM GSQRT
      CALL cdf_define(nwout, vn_bsupumnc, bsupumnc, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmnc, bsupvmnc, dimname=r3dim)
!     IF (lfreeb) THEN
!         CALL cdf_define(nwout, vn_rbc, rbc, 
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_rbc, ln_rbc, units='m')
!         CALL cdf_define(nwout, vn_zbs, zbs, 
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_zbs, ln_zbs, units='m')
!        IF (lasym) THEN
!           CALL cdf_define(nwout, vn_rbs, rbs, 
!    1                dimname=(/'n_mode','m_mode'/))
!           CALL cdf_define(nwout, vn_zbc, zbc, 
!    1                dimname=(/'n_mode','m_mode'/))
!        END IF
!     END IF

      IF (.NOT. lasym) GO TO 800

      ALLOCATE (gmns(mnmax_nyq,ns), bmns(mnmax_nyq,ns),
     1          bsubumns(mnmax_nyq,ns), bsubvmns(mnmax_nyq,ns),
     2          bsubsmnc(mnmax_nyq,ns), bsupumns(mnmax_nyq,ns),
     3          bsupvmns(mnmax_nyq,ns))
      CALL cdf_define(nwout, vn_rmns, rmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_rmns, ln_rmns, units='m')
      CALL cdf_define(nwout, vn_zmnc, zmnc, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_zmnc, ln_zmnc, units='m')
      CALL cdf_define(nwout, vn_lmnc, lmnc, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_lmnc, ln_lmnc)
      CALL cdf_define(nwout, vn_gmns, gmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_gmns, ln_gmns)
      CALL cdf_define(nwout, vn_bmns, bmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bmns, ln_bmns)
      CALL cdf_define(nwout, vn_bsubumns, bsubumns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubumns, ln_bsubumns)
      CALL cdf_define(nwout, vn_bsubvmns, bsubvmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubvmns, ln_bsubvmns)
      CALL cdf_define(nwout, vn_bsubsmnc, bsubsmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubsmnc, ln_bsubsmnc)

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_define(nwout, vn_bsupumns, bsupumns, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmns, bsupvmns, dimname=r3dim)

 800  CONTINUE

!================================
! Write Variables
!================================

! Scalars
      CALL cdf_write(nwout, vn_version, vversion)
      CALL cdf_write(nwout, vn_extension, input_extension)
      CALL cdf_write(nwout, vn_mgrid, mgrid_file)
      CALL cdf_write(nwout, vn_magen, wb)
      CALL cdf_write(nwout, vn_therm, wp)
      CALL cdf_write(nwout, vn_gam, gamma)
      CALL cdf_write(nwout, vn_maxr, rmax_surf)
      CALL cdf_write(nwout, vn_minr, rmin_surf)
      CALL cdf_write(nwout, vn_maxz, zmax_surf)
      CALL cdf_write(nwout, vn_fp, nfp)
      CALL cdf_write(nwout, vn_radnod, ns)
      CALL cdf_write(nwout, vn_polmod, mpol)
      CALL cdf_write(nwout, vn_tormod, ntor)
      CALL cdf_write(nwout, vn_maxmod, mnmax)
      CALL cdf_write(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_write(nwout, vn_maxit, niter)
      CALL cdf_write(nwout, vn_actit, itfsq)
      CALL cdf_write(nwout, vn_asym, lasym)
      CALL cdf_write(nwout, vn_recon, lrecon)
      CALL cdf_write(nwout, vn_free, lfreeb)
      CALL cdf_write(nwout, vn_error, ier_flag)
      CALL cdf_write(nwout, vn_aspect, aspect)
      CALL cdf_write(nwout, vn_beta, betatot)
      CALL cdf_write(nwout, vn_pbeta, betapol)
      CALL cdf_write(nwout, vn_tbeta, betator)
      CALL cdf_write(nwout, vn_abeta, betaxis)
      CALL cdf_write(nwout, vn_b0, b0)
      CALL cdf_write(nwout, vn_rbt0, rbtor0)
      CALL cdf_write(nwout, vn_rbt1, rbtor)
      CALL cdf_write(nwout, vn_sgs, NINT(signgs))
      CALL cdf_write(nwout, vn_lar, IonLarmor)
      CALL cdf_write(nwout, vn_modB, volAvgB)
      CALL cdf_write(nwout, vn_ctor, ctor/mu0)
      CALL cdf_write(nwout, vn_amin, Aminor_p)
      CALL cdf_write(nwout, vn_rmaj, Rmajor_p)
      CALL cdf_write(nwout, vn_vol, volume_p)

      IF (lrecon) THEN
         CALL cdf_write(nwout, vn_mse, imse2-1)
         CALL cdf_write(nwout, vn_thom, itse)
      END IF

      CALL cdf_write(nwout, vn_nextcur, nextcur)
      IF (nextcur .gt. 0) THEN
         CALL cdf_write(nwout, vn_extcur, extcur(1:nextcur))
         CALL cdf_write(nwout, vn_mgmode, mgrid_mode)
      ENDIF
      IF (lfreeb) THEN
         CALL cdf_write(nwout, vn_flp, nobser)
         CALL cdf_write(nwout, vn_nobd, nobd)
         CALL cdf_write(nwout, vn_nbset, nbsets)
         IF (nextcur .gt. 0)
     1   CALL cdf_write(nwout, vn_curlab, curlabel(1:nextcur))
      END IF

! 1D Arrays
      IF (nbsets .gt. 0) CALL cdf_write(nwout,vn_nbfld,nbfld(1:nbsets))

      IF (.not.lwrite) GO TO 1000

      CALL cdf_write(nwout, vn_pmod, xm)
      CALL cdf_write(nwout, vn_tmod, xn)
      CALL cdf_write(nwout, vn_pmod_nyq, xm_nyq)
      CALL cdf_write(nwout, vn_tmod_nyq, xn_nyq)

      RADIUS1: DO js = 1, ns

         CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1, xc, js)

         rmnc(:,js) = rmnc1(:)
         zmns(:,js) = zmns1(:)
         lmns(:,js) = lmns1(:)
         IF (lasym) THEN
            rmns(:,js) = rmns1(:)
            zmnc(:,js) = zmnc1(:)
            lmnc(:,js) = lmnc1(:)
         END IF

      END DO RADIUS1

!
!     INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY
!
      WHERE (NINT(xm) .eq. 1) lmns(:,1) = lmns(:,2)
      DO js = ns,2,-1
         WHERE (MOD(NINT(xm),2) .eq. 0) 
            lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
         ELSE WHERE
            lmns(:,js) = p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
         END WHERE
      END DO

      lmns(:,1) = 0  
      raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
      zaxis_cs(0:ntor) = zmns(1:ntor+1,1)

      CALL cdf_write(nwout, vn_racc, raxis_cc(0:ntor))
      CALL cdf_write(nwout, vn_zacs, zaxis_cs(0:ntor)) 
      CALL cdf_write(nwout, vn_rmnc, rmnc)
      CALL cdf_write(nwout, vn_zmns, zmns)
      CALL cdf_write(nwout, vn_lmns, lmns)

      IF (.not.lasym) GOTO 900

      WHERE (NINT(xm) .eq. 1) lmns(:,1) = lmns(:,2)
      DO js = ns,2,-1
         WHERE (MOD(NINT(xm),2) .eq. 0) 
            lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
         ELSE WHERE
            lmnc(:,js) = p5*(sm(js)*lmnc(:,js) + sp(js-1)*lmnc(:,js-1))
         END WHERE
      END DO

      lmnc(:,1) = 0;   
      raxis_cs(0:ntor) = rmns(1:ntor+1,1)
      zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)
      CALL cdf_write(nwout, vn_racs, raxis_cs(0:ntor))
      CALL cdf_write(nwout, vn_zacc, zaxis_cc(0:ntor)) 
      CALL cdf_write(nwout, vn_rmns, rmns)
      CALL cdf_write(nwout, vn_zmnc, zmnc)
      CALL cdf_write(nwout, vn_lmnc, lmnc)

 900  CONTINUE

!     NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
      IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
      IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)
      tmult = p5/r0scale**2
      IF (lasym) tmult = tmult*p5

      RADIUS2: DO js = 2, ns
         bmod = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
         gmn = 0
         bmn = 0
         bsubumn = 0
         bsubvmn = 0
         bsubsmn = 0
         bsupumn = 0
         bsupvmn = 0
         MN2: DO mn = 1, mnmax_nyq
            n = NINT(xn_nyq(mn))/nfp
            m = NINT(xm_nyq(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1 
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))          !cos(mu - nv)
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))          !sin(mu - nv)
                  bmn(mn) = bmn(mn) + tcosi*bmod(lk)
                  gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lk)
                  bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
                  bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
                  bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
                  bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk)
                  bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk)
               END DO
            END DO
         END DO MN2

!        Do u -> -u integral, 0, pi ; flip sign on sin theta terms
         IF (lasym) THEN
         MN2_ASYM: DO mn = 1, mnmax_nyq
            n = NINT(xn_nyq(mn))/nfp
            m = NINT(xm_nyq(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1 
                  lr = uminus(lk)
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) -
     1                       sgn*sinmui(j,m)*sinnv(k,n1))         !cos(-mu-nv)
                  tsini =-dmult*(sinmui(j,m)*cosnv(k,n1) +
     1                       sgn*cosmui(j,m)*sinnv(k,n1))         !sin(-mu-nv)
                  bmn(mn) = bmn(mn) + tcosi*bmod(lr)
                  gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lr)
                  bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lr)
                  bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lr)
                  bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lr)
                  bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lr)
                  bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lr)
               END DO
            END DO
         END DO MN2_ASYM
         END IF
         
         IF (js .eq. ns/2) bmodmn = bmn
         IF (js .eq. ns) bmodmn1 = bmn
         gmnc(:,js) = gmn(:)
         bmnc(:,js) = bmn(:)
         bsubumnc(:,js) = bsubumn(:)
         bsubvmnc(:,js) = bsubvmn(:)
         bsubsmns(:,js) = bsubsmn(:)
         bsupumnc(:,js) = bsupumn(:)
         bsupvmnc(:,js) = bsupvmn(:)

      END DO RADIUS2

      gmnc(:,1) = 0; bmnc(:,1) = 0; bsubumnc(:,1) = 0
      bsubvmnc(:,1) = 0
      bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3)
      bsupumnc(:,1) = 0;  bsupvmnc(:,1) = 0
      CALL cdf_write(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubsmns, bsubsmns)      !Full mesh
!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_write(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_write(nwout, vn_bsupvmnc, bsupvmnc)

      DEALLOCATE (gmnc, bmnc, bsubumnc, bsubvmnc,
     2            bsubsmns, bsupumnc, bsupvmnc)

      IF (.not.lasym) GO TO 200

      RADIUS3: DO js = 2, ns
         bmod = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
         gmn = 0
         bmn = 0
         bsubumn = 0
         bsubvmn = 0
         bsubsmn = 0
         bsupumn = 0
         bsupvmn = 0
         MN3: DO mn = 1, mnmax_nyq
            n = NINT(xn_nyq(mn))/nfp
            m = NINT(xm_nyq(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))
                  bmn(mn) = bmn(mn) + tsini*bmod(lk)
                  gmn(mn) = gmn(mn) + tsini*gsqrt(js,lk)
                  bsubumn(mn) = bsubumn(mn) + tsini*bsubu(js,lk)
                  bsubvmn(mn) = bsubvmn(mn) + tsini*bsubv(js,lk)
                  bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubs(js,lk)
                  bsupumn(mn) = bsupumn(mn) + tsini*bsupu(js,lk)
                  bsupvmn(mn) = bsupvmn(mn) + tsini*bsupv(js,lk)
               END DO
            END DO
         END DO MN3
         MN3_ASYM: DO mn = 1, mnmax_nyq
            n = NINT(xn_nyq(mn))/nfp
            m = NINT(xm_nyq(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  lr = uminus(lk)
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) -
     1                       sgn*sinmui(j,m)*sinnv(k,n1))
                  tsini =-dmult*(sinmui(j,m)*cosnv(k,n1) +
     1                       sgn*cosmui(j,m)*sinnv(k,n1))
                  bmn(mn) = bmn(mn) + tsini*bmod(lr)
                  gmn(mn) = gmn(mn) + tsini*gsqrt(js,lr)
                  bsubumn(mn) = bsubumn(mn) + tsini*bsubu(js,lr)
                  bsubvmn(mn) = bsubvmn(mn) + tsini*bsubv(js,lr)
                  bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubs(js,lr)
                  bsupumn(mn) = bsupumn(mn) + tsini*bsupu(js,lr)
                  bsupvmn(mn) = bsupvmn(mn) + tsini*bsupv(js,lr)
               END DO
            END DO
         END DO MN3_ASYM
   
         gmns(:,js) = gmn(:)
         bmns(:,js) = bmn(:)
         bsubumns(:,js) = bsubumn(:)
         bsubvmns(:,js) = bsubvmn(:)
         bsubsmnc(:,js) = bsubsmn(:)
         bsupumns(:,js) = bsupumn(:)
         bsupvmns(:,js) = bsupvmn(:)

      END DO RADIUS3

      gmns(:,1) = 0; bmns(:,1) = 0; bsubumns(:,1) = 0
      bsubvmns(:,1) = 0
      bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
      bsupumns(:,1) = 0;  bsupvmns(:,1) = 0
      CALL cdf_write(nwout, vn_gmns, gmns)
      CALL cdf_write(nwout, vn_bmns, bmns) 
      CALL cdf_write(nwout, vn_bsubumns, bsubumns)
      CALL cdf_write(nwout, vn_bsubvmns, bsubvmns)
      CALL cdf_write(nwout, vn_bsubsmnc, bsubsmnc)
!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_write(nwout, vn_bsupumns, bsupumns)
      CALL cdf_write(nwout, vn_bsupvmns, bsupvmns)

      DEALLOCATE (gmns, bmns, bsubumns, bsubvmns,
     2            bsubsmnc, bsupumns, bsupvmns)

 200  CONTINUE

      IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
      IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)

!     FULL-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
      CALL cdf_write(nwout, vn_iotaf, iotaf(1:ns)) 
      CALL cdf_write(nwout, vn_presf, presf/mu0) 
      CALL cdf_write(nwout, vn_phi, phi) 
      CALL cdf_write(nwout, vn_phipf, phipf)
      CALL cdf_write(nwout, vn_jcuru, jcuru/mu0)
      CALL cdf_write(nwout, vn_jcurv, jcurv/mu0)
      CALL cdf_write(nwout, vn_jdotb, jdotb)
      CALL cdf_write(nwout, vn_bgrv, bdotgradv)
 
!     HALF-MESH quantities
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0; 
      buco(1) = 0; bvco(1) = 0; vp(1) = 0; overr(1) = 0;  specw(1) = 1
      beta_vol(1) = 0
      CALL cdf_write(nwout, vn_iotah, iotas(1:ns))
      CALL cdf_write(nwout, vn_mass, mass/mu0) 
      CALL cdf_write(nwout, vn_presh, pres(1:ns)/mu0)
      CALL cdf_write(nwout, vn_betah, beta_vol)
      CALL cdf_write(nwout, vn_buco, buco)
      CALL cdf_write(nwout, vn_bvco, bvco) 
      CALL cdf_write(nwout, vn_vp, vp(1:ns))
      CALL cdf_write(nwout, vn_specw, specw)
      CALL cdf_write(nwout, vn_phip, phips(1:ns))
      CALL cdf_write(nwout, vn_overr, overr(1:ns))

!     MERCIER_CRITERION
      CALL cdf_write(nwout, vn_merc, Dmerc)
      CALL cdf_write(nwout, vn_mshear, Dshear)
      CALL cdf_write(nwout, vn_mwell, Dwell)
      CALL cdf_write(nwout, vn_mcurr, Dcurr)
      CALL cdf_write(nwout, vn_mgeo, Dgeod)
      CALL cdf_write(nwout, vn_equif, equif)

      CALL cdf_write(nwout, vn_fsq, fsqt(1:nstore_seq))
      CALL cdf_write(nwout, vn_wdot, wdot(1:nstore_seq))  

!-----------------------------------------------
!     DATA AND MSE FITS : HAVE NOT CONVERTED TO NETCDF
!     SINCE THIS WILL BE REPLACED SOON
!-----------------------------------------------
      IF (.not.lrecon) GOTO 950

 950  CONTINUE
!-----------------------------------------------
!     FREE BOUNDARY DATA
!-----------------------------------------------
      CALL freeb_data(rmnc1, zmns1, rmns1, zmnc1, bmodmn, bmodmn1)

 1000 CONTINUE

      rzl_array = 0

      CALL cdf_close(nwout)

      END SUBROUTINE wrout_cdf
