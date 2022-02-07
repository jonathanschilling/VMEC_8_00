      MODULE read_boozer_mod
!
!     USAGE:
!
!     Use READ_BOOZ_MOD to include variables dynamically ALlocated
!     in the module
!     Call DEALLOCATE_READ_BOOZER to free this memory when it is no longer needed
!
      USE stel_kinds
      IMPLICIT NONE
!DEC$ IF DEFINED (NETCDF)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER ::
     1 vn_nfp="nfp_b", vn_ns="ns_b", vn_aspect="aspect_b",
     2 vn_rmax="rmax_b", vn_rmin="rmin_b", vn_betaxis="betaxis_b",
     3 vn_mboz="mboz_b", vn_nboz="nboz_b", vn_mnboz="mnboz_b",
     4 vn_version="version", vn_iota="iota_b", vn_pres="pres_b",
     5 vn_beta="beta_b", vn_phip="phip_b", vn_phi="phi_b",
     6 vn_bvco="bvco_b", vn_buco="buco_b", vn_ixm="ixm_b",
     7 vn_ixn="ixn_b", vn_bmn="bmn_b", vn_rmnc="rmnc_b",
     8 vn_zmns="zmns_b", vn_pmns="pmns_b", vn_gmn="gmn_b",
     9 vn_jlist="jlist"

!DEC$ ENDIF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mnboz_b, mboz_b, nboz_b, nfp_b, ns_b
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx_b, ixm_b, ixn_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: iota_b, pres_b,
     1    phip_b, phi_b, beta_b, buco_b, bvco_b
      REAL(rprec) :: aspect_b, rmax_b, rmin_b, betaxis_b
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   bmn_b, rmnc_b, zmns_b, pmns_b, gmn_b, packed2d

      CONTAINS

      SUBROUTINE read_boozer_file (file_or_extension, ierr, iopen)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER*(*) :: file_or_extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_booz = 14
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
      LOGICAL :: isnc
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE BOOZMN FILE CREATED BY THE BOOZ_XFORM CODE
!     AND STORES THE DATA IN THE READ_BOOZ_MOD MODULE
!
!     CHECK FOR netcdf FILE EXTENSION (*.nc)
!
      filename = 'boozmn'
      CALL parse_extension(filename, file_or_extension, isnc)

      IF (isnc) THEN
!DEC$ IF DEFINED (NETCDF)
         CALL read_boozer_nc(filename, ierr)
!DEC$ ELSE
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
!DEC$ ENDIF
      ELSE
         iunit = unit_booz
         CALL safe_open (iunit, ierr, filename, 'old', 'unformatted')
         IF (ierr .eq. 0) CALL read_boozer_bin(iunit, ierr)
         CLOSE(unit=iunit)
      END IF

      IF (PRESENT(iopen)) iopen = ierr

      END SUBROUTINE read_boozer_file

!DEC$ IF DEFINED (NETCDF)
      SUBROUTINE read_boozer_nc(filename, ierr)
      USE stel_constants, ONLY: zero
      USE ezcdf
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      CHARACTER*(*) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER :: nbooz, nsval, ilist
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlist
      CHARACTER*38 :: version
C-----------------------------------------------
! Open cdf File
      call cdf_open(nbooz,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening boozmn .nc file'
         RETURN
      END IF

! Read in scalar variables
      CALL cdf_read(nbooz, vn_nfp, nfp_b)
      CALL cdf_read(nbooz, vn_ns, ns_b)
      CALL cdf_read(nbooz, vn_aspect, aspect_b)
      CALL cdf_read(nbooz, vn_rmax, rmax_b)
      CALL cdf_read(nbooz, vn_rmin, rmin_b)
      CALL cdf_read(nbooz, vn_betaxis, betaxis_b)
      CALL cdf_read(nbooz, vn_mboz, mboz_b)
      CALL cdf_read(nbooz, vn_nboz, nboz_b)
      CALL cdf_read(nbooz, vn_mnboz, mnboz_b)
      CALL cdf_read(nbooz, vn_version, version)

!  1D arrays (skip inquiry statements for now, assume correct in file)
      IF (ALLOCATED(iota_b)) CALL read_boozer_deallocate
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),
     1  phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), 
     2  ixm_b(mnboz_b), ixn_b(mnboz_b), stat=ierr)
      IF (ierr .ne. 0) THEN
        PRINT *,' Allocation error in read_boozer_file'
        RETURN
      END IF

      CALL cdf_read(nbooz, vn_iota, iota_b)
      CALL cdf_read(nbooz, vn_pres, pres_b)
      CALL cdf_read(nbooz, vn_beta, beta_b)
      CALL cdf_read(nbooz, vn_phip, phip_b)
      CALL cdf_read(nbooz, vn_phi,  phi_b)
      CALL cdf_read(nbooz, vn_bvco, bvco_b)
      CALL cdf_read(nbooz, vn_buco, buco_b)
      CALL cdf_read(nbooz, vn_ixm, ixm_b)
      CALL cdf_read(nbooz, vn_ixn, ixn_b)

      CALL cdf_inquire(nbooz, vn_jlist, dimlens)
      ALLOCATE (jlist(1:dimlens(1)), stat=ierr)
      CALL cdf_read(nbooz, vn_jlist, jlist)

      idx_b = 0
      ilist = SIZE(jlist)
      DO ilist = 1, SIZE(jlist)
         nsval = jlist(ilist)
         idx_b(nsval) = 1
      END DO

!  2D arrays
      ALLOCATE (bmn_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),
     1  zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b), gmn_b(mnboz_b,ns_b),
     2  packed2d(mnboz_b, ilist), stat = ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

!
!  Note: Must unpack these 2D arrays, only jlist-ed radial nodes store in file
!
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmn_b = 0; gmn_b = 0
      CALL unpack_cdf(nbooz, vn_bmn, bmn_b)
      CALL unpack_cdf(nbooz, vn_rmnc, rmnc_b)
      CALL unpack_cdf(nbooz, vn_zmns, zmns_b)
      CALL unpack_cdf(nbooz, vn_pmns, pmns_b)
      CALL unpack_cdf(nbooz, vn_gmn, gmn_b)

      DEALLOCATE (jlist, packed2d)

! Close cdf File
      CALL cdf_close(nbooz, ierr)

      END SUBROUTINE read_boozer_nc
!DEC$ ENDIF

      SUBROUTINE read_boozer_bin(iunit, ierr)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr, iunit
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsval, jsize, js
      CHARACTER*38 :: version
C-----------------------------------------------

      READ(iunit, iostat=ierr, err=100) nfp_b, ns_b, aspect_b,
     1   rmax_b, rmin_b, betaxis_b

      IF (ALLOCATED(iota_b)) CALL read_boozer_deallocate
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),
     1  phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), stat=ierr)
      IF (ierr .ne. 0) THEN
        PRINT *,' Allocation error in read_boozer_file'
        RETURN
      END IF
 
      iota_b(1) = 0; pres_b(1) = 0; beta_b(1) = 0
      phip_b(1) = 0; phi_b(1) = 0; bvco_b(1) = 0
      buco_b(1) = 0

      DO nsval = 2, ns_b
         READ(iunit, iostat=ierr, err=100) iota_b(nsval),
     1   pres_b(nsval), beta_b(nsval), phip_b(nsval), phi_b(nsval),
     2   bvco_b(nsval), buco_b(nsval)
      END DO

      READ(iunit, iostat=ierr, err=100) mboz_b, nboz_b, mnboz_b, jsize
      READ(iunit, iostat=js) version

      ALLOCATE (bmn_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),
     1  zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b), gmn_b(mnboz_b,ns_b),
     2  ixm_b(mnboz_b), ixn_b(mnboz_b), stat = ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

!
!     idx_b:  = 0, surface data not in file; = 1, surface data in file
!
      idx_b = 0
      ixm_b = 0
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmn_b = 0; gmn_b = 0

      READ(iunit,iostat=ierr,err=100) ixn_b(:mnboz_b), ixm_b(:mnboz_b)

      DO js = 1, jsize
        READ(iunit, iostat=ierr, END=200, err=100) nsval
        IF ((nsval.gt.ns_b) .or. (nsval.le.0)) CYCLE

        idx_b(nsval) = 1

        READ(iunit,iostat=ierr,err=100, END=200) bmn_b(:mnboz_b,nsval),
     1       rmnc_b(:mnboz_b,nsval), zmns_b(:mnboz_b,nsval), 
     2       pmns_b(:mnboz_b,nsval), gmn_b(:mnboz_b,nsval)
      END DO

 100  CONTINUE
      IF (ierr .gt. 0)
     1    PRINT *,' Error reading in subroutine read_boozer_file:',
     2            ' ierr = ', ierr
 200  CONTINUE
      IF (ierr .lt. 0) ierr = 0       !End-of-file, ok
      CLOSE(iunit)

      END SUBROUTINE read_boozer_bin


      SUBROUTINE read_boozer_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------

      IF (ALLOCATED(iota_b)) DEALLOCATE (iota_b, pres_b, beta_b,
     1    phip_b, phi_b, bvco_b, buco_b, idx_b, stat = istat)

      IF (ALLOCATED(bmn_b)) DEALLOCATE (bmn_b, rmnc_b,
     1   zmns_b, pmns_b, gmn_b, ixm_b, ixn_b, stat = istat)

      END SUBROUTINE read_boozer_deallocate


!DEC$ IF DEFINED (NETCDF)
      SUBROUTINE unpack_cdf (nbooz, var_name, array2d)
      USE stel_kinds, ONLY: rprec
      USE ezcdf
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nbooz
      INTEGER :: nsval, icount
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: array2d
      CHARACTER*(*), INTENT(in) :: var_name
C-----------------------------------------------
!
!     Read into temporary packed array, packed2d
!
      CALL cdf_read(nbooz, var_name, packed2d)

      array2d = 0; icount = 1

      DO nsval = 1, ns_b
         IF (idx_b(nsval) .eq. 1) THEN
            array2d(:,nsval) = packed2d(:,icount)
            icount = icount + 1
         END IF
      END DO

      END SUBROUTINE unpack_cdf
!DEC$ ENDIF


      END MODULE read_boozer_mod


 