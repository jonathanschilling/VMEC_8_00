      SUBROUTINE vsetup (lreset, iseq_count)
      USE vmec_main
      USE vacmod
      USE realspace
      USE vsvd
      USE timer_sub
      USE mgrid_mod, ONLY: nbcoil_max, nlim_max, nextcur, mgrid_mode
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iseq_count
      LOGICAL :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: mpol_default = 6
      INTEGER, PARAMETER :: ntor_default = 0
C-----------------------------------------------
!
!     Default and initial values
!
      lasym = .false.
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0
      nzeta  = 0
      nthreed = 0
      lfreeb = .true.              !reset later after input file READ
      lrecon = .true.
      loldout = .false.
      laddout = .false.
      ldiagno = .false.
      lmoreiter = .true.
      lmac = .false.
      ledge_dump = .false.
      liota = .false.

      ipmass = 0     ! profile selection 0=original polynomial
      ipiota = 0     ! profile selection 0=original polynomial
      ipcurr = 0     ! profile selection 0=original polynomial

      z00 = zero
      mgrid_file = 'NONE'
      nextcur = 0
      extcur = 0
      mgrid_mode = 'S'             !Assume scaled mode; read in from mgrid in free-bdy mode

!
!     ZERO ARRAYS WHICH MAY BE READ IN
!
      rbc = zero
      rbs = zero
      zbc = zero
      zbs = zero

      am = zero
      ai = zero
      ac = cbig
      aphi = zero
      extcur = zero
      bloat = one
      pres_scale = one

      ns_array = 0
      ftol_array = zero
      raxis_cc = zero;  raxis_cs = zero
      zaxis_cc = zero;  zaxis_cs = zero
      gamma = zero
      spres_ped = one

      iequi = 0
      ivac  = -1
      delbsq = one
      delt = 1.1
      tcon0 = one
      curtor = 1.e30_dp
      time_slice = zero
      IF (lreset) THEN
         pfac   = one
         phifac = one
         timer = zero
      ENDIF
      fsqr = one
      fsqz = one
      ftolv = fsqr
!
!     Reconstruction stuff
!
      lpofr = .true.
      lpprof = .true.
      icurrout = 0
      total_chi_square_n = one
      total_chisq_n0 = total_chi_square_n
      nchistp = 0
      imse = -1
      itse = 0
      isnodes = 0
      ipnodes = 0
      iopt_raxis = 1
      imatch_phiedge = 1
      nflxs = 0
      nbfld = 0
      mseangle_offset = zero
      mseangle_offsetm = zero
      pres_offset = zero
      sigma_current = 1.e30_dp
      sigma_delphid = 1.e30_dp
      tensi = one
      tensp = one
      tensi2 = zero
      fpolyi = one
      presfac = one
      phidiam = 1.e30_dp

      mseprof = one
      indxflx = 0
      indxbfld = 0
      sigma_stark = 1.1*cbig
      sigma_thom = 1.1*cbig
      sigma_flux = 1.1*cbig
      sigma_b = 1.1*cbig
!
!     FREE-BOUNDARY STUFF, READ IN FIRST TIME ONLY
!
      IF (iseq_count .eq. 0) THEN
        nbcoil_max = 0
        nlim_max = 0
      END IF

      END SUBROUTINE vsetup
