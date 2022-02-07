      MODULE vmec_params
      USE stel_kinds, ONLY: rprec1 => rprec, dp1 => dp
      USE vparams, ONLY: mpold1 => mpold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: meven = 0, modd = 1
      INTEGER, PARAMETER :: ndamp = 10
      INTEGER, PARAMETER :: ns4 = 25

      INTEGER, PRIVATE :: ink
      INTEGER, PARAMETER, DIMENSION(0:mpold1) ::
     1  jmin1 = (/ 1,1,(2,ink=2,mpold1) /),        !starting js(m) values where R,Z are non-zero
     2  jmin2 = (/ 1,2,(2,ink=2,mpold1) /),        !starting js(m) values for which R,Z are evolved
!     3  jlam  = (/ 1,2,(2,ink=2,mpold1) /)         !starting js(m) values for which Lambda is evolved
     3  jlam  = (/ 2,2,(2,ink=2,mpold1) /)         !starting js(m) values for which Lambda is evolved

      INTEGER, PARAMETER :: norm_term_flag=0, bad_init_jac_flag=1, 
     1                      more_iter_flag=4, bad_jacobian_flag=6
     
      REAL(rprec1), PARAMETER :: pdamp = 0.05
      CHARACTER*(*), PARAMETER :: version_ = '8.00'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ntmax, rcc, rss, rsc, rcs, zsc, zcs, zcc, zss
      INTEGER :: mnyq, nnyq
      INTEGER, ALLOCATABLE :: uminus(:)
      REAL(rprec1), ALLOCATABLE :: mscale(:), nscale(:)
      REAL(rprec1) :: signgs
!-----------------------------------------------

      END MODULE vmec_params
