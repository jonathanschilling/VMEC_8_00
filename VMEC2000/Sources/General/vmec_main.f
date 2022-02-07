      MODULE vmec_main
      USE vmec_dim
      USE vmec_input
      USE vmec_persistent
      USE vmec_params, ONLY: ndamp
      USE vparams
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!     lthreed:  .true., 3d plasma; .false., 2d plasma (n=0)
!     liota:    .true., iotaf determined from Force-iota ~ Bsubu - icurv
!               .false.,iotas solved for and iotaf obtained by interpolation (default)
!
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    ard, arm, brd, brm, azd, azm, bzd, bzm, bmin, bmax
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    crd, iotaf, phipf, mass, phi, presf, beta_vol, jcuru, jcurv, 
     2    jdotb, buco, bvco, bdotgradv, equif, specw, tcon, 
     3    psi, yellip, yinden, ytrian, yshift, ygeo, overr, 
     4    sm, sp, iotas, phips, pres, vp, jpar2, jperp2, bdotb, 
     5    blam, clam, dlam, icurv, vpphi, presgrad,
     6    r01, z01
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam
      REAL(rprec), DIMENSION(0:mpol1d,3) :: xmpq
      REAL(rprec), DIMENSION(0:mpol1d) :: faccon
      REAL(rprec) :: dcon, currv, aspect, hs, ohs, voli, 
     1   signiota, rc0mse, r00, r0scale, z00, dkappa, fsqsum0,
     2   pressum0, fnorm, fsqr, fsqz, fsql, fnorm1, fsqr1, fsqz1, 
     3   fsql1, fsq, fedge, wb, wp, r00b, z00b, fz00_edge
      REAL(rprec), DIMENSION(nstore_seq) :: fsqt, wdot
      REAL(rprec) :: ftolv, otav
      REAL(rprec), DIMENSION(ndamp) :: otau
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::
     1    rmn_bdy, zmn_bdy
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsqsav
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubu0, dbsq, rbsq
      REAL(rprec) :: rbtor, rbtor0, ctor, delbsq
      REAL(rprec), DIMENSION(ndatafmax) ::
     1  spfa, spfa2, hp, sifa, sifa2, hi
      LOGICAL :: lthreed, liota
      INTEGER, DIMENSION(:), ALLOCATABLE :: ireflect
      INTEGER :: multi_ns_grid, iequi, irst,
     1    iter1, iter2, ijacob, itfsq, iresidue, neqs, neqs1,
     2    neqs2, irzloff, ivac, ndatap, ndatai
C-----------------------------------------------
      END MODULE vmec_main
