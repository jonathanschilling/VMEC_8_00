      MODULE bootsj_input
      USE stel_kinds
      INTEGER :: nrho, mbuse, nbuse, isymm0
      REAL(rprec) :: tempres, ate(0:11), ati(0:11)
      REAL(rprec) :: zeff1, dens0, teti, damp, damp_bs

      NAMELIST /bootin/ nrho, mbuse, nbuse, zeff1, dens0, teti, tempres,
     1   damp, damp_bs, isymm0, ate, ati

      CONTAINS

      SUBROUTINE read_boot_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=bootin, iostat=istat)

      END SUBROUTINE read_boot_namelist

      END MODULE bootsj_input
