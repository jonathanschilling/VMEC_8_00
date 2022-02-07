      MODULE fdjac_mod
      USE stel_kinds

      INTEGER, PARAMETER :: flag_singletask = -1, flag_cleanup = -100

      INTEGER :: m, n, ncnt, max_processors, num_lm_params
      REAL(rprec), DIMENSION(:), POINTER :: xp, wap
      LOGICAL, DIMENSION(:), ALLOCATABLE :: flip
      REAL(rprec) :: eps

      END MODULE fdjac_mod
