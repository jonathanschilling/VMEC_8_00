      MODULE vmec_input
      USE vparams, ONLY: rprec, dp, mpol1d, ntord, ndatafmax
      USE vsvd0
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      dsiobt:   measured flux loop signals corresponding to the
!                combination of signals in iconnect array
!     indxflx:   array giving index of flux measurement in iconnect array
!    indxbfld:   array giving index of bfield measurement used in matching
!    bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
!                the orientation br*cos(abcoil) + bz*sin(abcoil)
!                of the m-th coil in the n-th set from mgrid file.
!       nflxs:   number of flux loop measurements used in matching
!    nbfld(n):   number of selected external bfield measurements in set n from nml file
!    pres_scale: factor used to scale pressure profile (default value = 1)
!                useful so user can fix profile and change beta without having to change
!                all AM coefficients separately
      INTEGER :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor,
     1           ntheta, nzeta, ipmass, ipiota, ipcurr
      INTEGER, DIMENSION(100) :: ns_array
      INTEGER :: imse, isnodes, itse, ipnodes, iopt_raxis,
     1   imatch_phiedge, nflxs
      INTEGER, DIMENSION(nbsetsp) :: nbfld
      INTEGER, DIMENSION(nfloops) :: indxflx
      INTEGER, DIMENSION(nbcoilsp,nbsetsp) :: indxbfld
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbs, zbc, rbc, zbs
      REAL(rprec) :: time_slice, curtor, delt, ftol, tcon0,
     1   gamma, phiedge, phidiam, sigma_current, sigma_delphid, tensi,
     2   tensp, tensi2, fpolyi, presfac, mseangle_offset, pres_offset,
     3   mseangle_offsetm, spres_ped, bloat, pres_scale
!     REAL(rprec), DIMENSION(0:10) :: am, ai, ac, aphi
      REAL(rprec), DIMENSION(0:12) :: am, ai, ac, aphi  !change by J.Geiger
      REAL(rprec), DIMENSION(0:ntord) :: raxis, zaxis                !!Backwards compatibility: Obsolete
      REAL(rprec), DIMENSION(0:ntord) :: raxis_cc, raxis_cs,
     1                                   zaxis_cc, zaxis_cs
      REAL(rprec), DIMENSION(100) :: ftol_array
      REAL(rprec), DIMENSION(nigroup) :: extcur
      REAL(rprec), DIMENSION(nmse) :: mseprof
      REAL(rprec), DIMENSION(ntse) :: rthom, datathom, sigma_thom
      REAL(rprec), DIMENSION(nmse) :: rstark, datastark,
     1    sigma_stark
      REAL(rprec), DIMENSION(nfloops) :: dsiobt, sigma_flux
      REAL(rprec), DIMENSION(nbcoilsp,nbsetsp) :: bbc, sigma_b
      REAL(rprec), DIMENSION(ndatafmax) :: psa, pfa, isa, ifa
      LOGICAL :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump,
     1           lasym, laddout, ldiagno, lmoreiter
     2         , lspectrum_dump, loptim           !!Obsolete
      CHARACTER*(100) :: mgrid_file
      CHARACTER*(120) :: arg1
      CHARACTER*(100) :: input_extension

      NAMELIST /indata/ mgrid_file, time_slice, nfp, ncurr, nsin,
     1   niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac, aphi,
     2   rbc, zbs, rbs, zbc, spres_ped, pres_scale, raxis_cc, zaxis_cs, 
     3   raxis_cs, zaxis_cc, mpol, ntor, ntheta, nzeta, 
     4   ns_array, ftol_array, tcon0, curtor, sigma_current, extcur,
     5   phiedge, psa, pfa, isa, ifa, imatch_phiedge, iopt_raxis, 
     6   tensi, tensp, mseangle_offset, mseangle_offsetm, imse, 
     7   isnodes, rstark, datastark, sigma_stark, itse, ipnodes, 
     8   presfac, pres_offset, rthom, datathom, sigma_thom, phidiam, 
     9   sigma_delphid, tensi2, fpolyi, nflxs, indxflx, dsiobt, 
     A   sigma_flux, nbfld, indxbfld, laddout, ldiagno, lmoreiter,
     A   bbc, sigma_b, lpofr, lfreeb, lrecon, lmac, loldout, lasym,
     B   ipmass, ipiota, ipcurr,
     B   ledge_dump, lspectrum_dump, loptim, bloat, raxis, zaxis

      NAMELIST /mseprofile/ mseprof

      CONTAINS

      SUBROUTINE read_indata_namelist (iunit, istat)
      REAL(rprec), PARAMETER :: smallno = -HUGE(1._dp)
      INTEGER :: iunit, istat

!
!     BACKWARDS COMPATIBILITY
!
      raxis = smallno
      zaxis = smallno
      
      READ (iunit, nml=indata, iostat=istat)

      WHERE (raxis .gt. smallno) raxis_cc = raxis
      WHERE (zaxis .gt. smallno) zaxis_cs = zaxis

      END SUBROUTINE read_indata_namelist

      SUBROUTINE read_mse_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=mseprofile, iostat=istat)

      END SUBROUTINE read_mse_namelist

      END MODULE vmec_input
