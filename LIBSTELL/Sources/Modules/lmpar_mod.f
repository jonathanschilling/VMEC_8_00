      MODULE lmpar_mod
      USE stel_kinds
      INTEGER :: nscan, ldfjac
      INTEGER, DIMENSION(:), POINTER :: ipvt
      REAL(rprec) :: pnorm, fnorm1, delta, par, spread_ratio
      REAL(rprec), DIMENSION(:), POINTER :: wa2p, wa3p, wa4p
      REAL(rprec), DIMENSION(:), POINTER :: diag, qtf
      REAL(rprec), DIMENSION(:,:), POINTER :: fjac
      LOGICAL :: lfirst_lm
      END MODULE lmpar_mod
