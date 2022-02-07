      MODULE fbal
      USE stel_kinds, ONLY: rprec1 => rprec
      LOGICAL, PARAMETER :: lforbal = .true.
      REAL(rprec1), DIMENSION(:), ALLOCATABLE :: rzu_fac, rru_fac,
     1  frcc_fac, fzsc_fac

      CONTAINS

      SUBROUTINE  calc_fbal(bsubu, bsubv)
      USE vmec_main, ONLY: buco, bvco, equif, 
     1                     jcurv, jcuru, iotaf, vp, pres, 
     2                     vpphi, presgrad, ohs
      USE vmec_params, ONLY: signgs
      USE vmec_dim, ONLY: ns, nrzt
      USE realspace, ONLY: wint, phip
C-----------------------------------------------
      REAL(rprec1), INTENT(in) :: bsubu(1:nrzt), bsubv(1:nrzt)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js
C-----------------------------------------------
      DO js = 2, ns
         buco(js) = SUM(bsubu(js:nrzt:ns)*wint(js:nrzt:ns))
         bvco(js) = SUM(bsubv(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO

      DO js = 2, ns-1
         jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
         jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))
         vpphi(js) = (vp(js+1)/phip(js+1) + vp(js)/phip(js))/2
         presgrad(js) = (pres(js+1) - pres(js))*ohs
         equif(js) = (-jcuru(js) + iotaf(js)*jcurv(js))/vpphi(js)
     1             + presgrad(js)
      END DO

      END SUBROUTINE calc_fbal

      END MODULE fbal
