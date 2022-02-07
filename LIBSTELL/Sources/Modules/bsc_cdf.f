!*******************************************************************************
!  File bsc_cdf.f
!  Contains the module bsc_cdf.f
!    Module for defining variables and writing netCDF files with the
!    derived types bsc_coil and bsc_coilcoll, from the bsc (Biot-Savart Coil) module
!
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module USEs the following modules:
!       stel_kinds
!       bsc
!       ezcdf
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  12.13.2002 - Initial Coding - Ed Lazarus, 
!  12.16.2002 - JDH Initial Comments, limit to bsc_cdf subroutines
!  12.17.2002 - JDH - return to using stel_kinds. Eliminate some unused variables.
!  12.18.2002 - JDH - Eliminated identifier. Added prefix. Made _coilcoll routines
!     call the _coil routines.
!  09.11.2003 - JDH - Added 'fil_rogo'wski information
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!*******************************************************************************

!*******************************************************************************
!  MODULE bsc_cdf
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      MODULE bsc_cdf

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------

      USE stel_kinds

!-------------------------------------------------------------------------------
!  Variable Names for netCDF
!-------------------------------------------------------------------------------

      CHARACTER (LEN=*), PARAMETER :: 
     &  vn_c_type = 'c_type',                                                  &
     &  vn_s_name = 's_name',                                                  &
     &  vn_l_name = 'l_name',                                                  &
     &  vn_current = 'current',                                                &         
     &  vn_raux = 'raux',                                                      &            
     &  vn_xnod = 'xnod',                                                      &            
     &  vn_ehnod = 'ehnod',                                                    &    
     &  vn_rcirc = 'rcirc',                                                    &        
     &  vn_xcent = 'xcent',                                                    &        
     &  vn_enhat = 'enhat',                                                    &
     &  vn_ave_n_area = 'ave_n_area'

!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINITION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_define_coil(this,lunit,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a bsc_coil

      USE bsc 
      USE ezcdf
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coil), INTENT (in) :: this
      INTEGER                      :: lunit
      CHARACTER (len=*)            :: prefix

!  this        bsc_coil - this is the coils that gets defined.
!  lunit       i/o unit number
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER (len=40) :: vn_c_type_use, vn_s_name_use,                      &
     &   vn_l_name_use, vn_current_use, vn_raux_use, vn_xnod_use,              &
     &   vn_rcirc_use, vn_xcent_use, vn_enhat_use, vn_ave_n_area_use
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      IF (LEN_TRIM(prefix) .GT. 0) THEN
         vn_c_type_use = bsc_cdf_mknam(prefix,vn_c_type)
         vn_s_name_use = bsc_cdf_mknam(prefix,vn_s_name)
         vn_l_name_use = bsc_cdf_mknam(prefix,vn_l_name)
         vn_current_use = bsc_cdf_mknam(prefix,vn_current)
         vn_raux_use = bsc_cdf_mknam(prefix,vn_raux)
         vn_xnod_use = bsc_cdf_mknam(prefix,vn_xnod)
         vn_rcirc_use = bsc_cdf_mknam(prefix,vn_rcirc)
         vn_xcent_use = bsc_cdf_mknam(prefix,vn_xcent)
         vn_enhat_use = bsc_cdf_mknam(prefix,vn_enhat)
         vn_ave_n_area_use = bsc_cdf_mknam(prefix,vn_ave_n_area)
      ELSE
         vn_c_type_use = vn_c_type
         vn_s_name_use = vn_s_name
         vn_l_name_use = vn_l_name
         vn_current_use = vn_current
         vn_raux_use = vn_raux
         vn_xnod_use = vn_xnod
         vn_rcirc_use = vn_rcirc
         vn_xcent_use = vn_xcent
         vn_enhat_use = vn_enhat
         vn_ave_n_area_use = vn_ave_n_area
      ENDIF      
         
! Define Components common to all c_types
      CALL cdf_define(lunit, TRIM(vn_c_type_use), this%c_type)
      CALL cdf_define(lunit, TRIM(vn_s_name_use), this%s_name)
      CALL cdf_define(lunit, TRIM(vn_l_name_use), this%l_name)
      CALL cdf_define(lunit, TRIM(vn_current_use), this%current)
      CALL cdf_define(lunit, TRIM(vn_raux_use), this%raux)

!  Particular coding, depending on c_type

      SELECT CASE (this%c_type)
      
      CASE ('fil_loop','floop') ! Filamentary Loop Variables
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_define(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED

      CASE ('fil_circ', 'fcirc') ! Filamentary Circle Variables
         CALL cdf_define(lunit, TRIM(vn_rcirc_use), this%rcirc)
         CALL cdf_define(lunit, TRIM(vn_xcent_use), this%xcent(1:3))
         CALL cdf_define(lunit, TRIM(vn_enhat_use), this%enhat(1:3))

      CASE ('fil_rogo') ! Rogowskis
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_define(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED
         CALL cdf_define(lunit, TRIM(vn_ave_n_area_use),                       &
     &         this%ave_n_area)
      
      END SELECT
      
      END SUBROUTINE bsc_cdf_define_coil

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_define_coilcoll(this,lunit)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coilcoll
!  To avoid duplicate names in the netCDF files, the variable names will
!  have a prefix added on.

      USE bsc 
      USE ezcdf
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coilcoll), INTENT (in) :: this
      INTEGER                          :: lunit

!  this        bsc_coilcoll - this is the coilcoll that gets defined.
!  lunit       i/o unit number
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER(iprec) :: i, n, ncoild
      INTEGER dimlens(2)
      CHARACTER*40 nowname
      CHARACTER (len=40) :: prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Not sure of the reason for the next IF test. JDH
!  
      IF (this%s_name .eq. ' ') THEN
         WRITE(*,*) 'this%s_name = one blank. bsc_cdf_define_coilcoll'
         WRITE(*,*) ' is returning'
         RETURN
      END IF
      
      ncoild = this%ncoil

!  Next loop could be augmented to make sure that the prefixes are unique.

      DO i = 1,ncoild  ! Loop over coils in the coilcoll
         prefix = this%coils(i)%s_name
         CALL bsc_cdf_define_coil(this%coils(i),lunit,prefix)
      END DO

      END SUBROUTINE bsc_cdf_define_coilcoll

!*******************************************************************************
! SECTION IV. WRITING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_write_coil(this,lunit,prefix)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coil

      USE bsc 
      USE ezcdf
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coil), INTENT (in) :: this
      INTEGER                      :: lunit
      CHARACTER (len=*)            :: prefix

!  this        bsc_coil - this is the coil that gets written.
!  lunit       i/o unit number
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER (len=40) :: vn_c_type_use, vn_s_name_use,                      &
     &   vn_l_name_use, vn_current_use, vn_raux_use, vn_xnod_use,              &
     &   vn_rcirc_use, vn_xcent_use, vn_enhat_use, vn_ave_n_area_use
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      IF (LEN_TRIM(prefix) .GT. 0) THEN
         vn_c_type_use = bsc_cdf_mknam(prefix,vn_c_type)
         vn_s_name_use = bsc_cdf_mknam(prefix,vn_s_name)
         vn_l_name_use = bsc_cdf_mknam(prefix,vn_l_name)
         vn_current_use = bsc_cdf_mknam(prefix,vn_current)
         vn_raux_use = bsc_cdf_mknam(prefix,vn_raux)
         vn_xnod_use = bsc_cdf_mknam(prefix,vn_xnod)
         vn_rcirc_use = bsc_cdf_mknam(prefix,vn_rcirc)
         vn_xcent_use = bsc_cdf_mknam(prefix,vn_xcent)
         vn_enhat_use = bsc_cdf_mknam(prefix,vn_enhat)
         vn_ave_n_area_use = bsc_cdf_mknam(prefix,vn_ave_n_area)
      ELSE
         vn_c_type_use = vn_c_type
         vn_s_name_use = vn_s_name
         vn_l_name_use = vn_l_name
         vn_current_use = vn_current
         vn_raux_use = vn_raux
         vn_xnod_use = vn_xnod
         vn_rcirc_use = vn_rcirc
         vn_xcent_use = vn_xcent
         vn_enhat_use = vn_enhat
         vn_ave_n_area_use = vn_ave_n_area
      ENDIF      
         
! Write Components common to all c_types
      CALL cdf_write(lunit, TRIM(vn_c_type_use), this%c_type)
      CALL cdf_write(lunit, TRIM(vn_s_name_use), this%s_name)
      CALL cdf_write(lunit, TRIM(vn_l_name_use), this%l_name)
      CALL cdf_write(lunit, TRIM(vn_current_use), this%current)
      CALL cdf_write(lunit, TRIM(vn_raux_use), this%raux)

!  Particular coding, depending on c_type

      SELECT CASE (this%c_type)
      
      CASE ('fil_loop','floop') ! Filamentary Loop Variables
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_write(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED

      CASE ('fil_circ', 'fcirc') ! Filamentary Circle Variables
         CALL cdf_write(lunit, TRIM(vn_rcirc_use), this%rcirc)
         CALL cdf_write(lunit, TRIM(vn_xcent_use), this%xcent(1:3))
         CALL cdf_write(lunit, TRIM(vn_enhat_use), this%enhat(1:3))

      CASE ('fil_rogo') ! Rogowskis
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_write(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED
         CALL cdf_write(lunit, TRIM(vn_ave_n_area_use),                        &
     &       this%ave_n_area)      
      
      END SELECT
      
      END SUBROUTINE bsc_cdf_write_coil

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_write_coilcoll(this,lunit)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coilcoll
!  To avoid duplicate names in the netCDF files, the variable names will
!  have a prefix added on.

      USE bsc 
      USE ezcdf
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coilcoll), INTENT (in) :: this
      INTEGER                          :: lunit

!  this        bsc_coilcoll - this is the coilcoll that gets written.
!  lunit       i/o unit number
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER(iprec) :: i, n, ncoild
      INTEGER dimlens(2)
      CHARACTER*40 nowname
      CHARACTER (len=40) :: prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Not sure of the reason for the next IF test. JDH
!  
      IF (this%s_name .eq. ' ') THEN
         WRITE(*,*) 'this%s_name = one blank. bsc_cdf_write_coilcoll'
         WRITE(*,*) ' is returning'
         RETURN
      END IF
      
      ncoild = this%ncoil

!  Next loop could be augmented to make sure that the prefixes are unique.

      DO i = 1,ncoild  ! Loop over coils in the coilcoll
         prefix = this%coils(i)%s_name
         CALL bsc_cdf_write_coil(this%coils(i),lunit,prefix)
      END DO

      END SUBROUTINE bsc_cdf_write_coilcoll

!*******************************************************************************
! SECTION IV. AUXILIARY FUNCTIONS
!*******************************************************************************
      
      FUNCTION bsc_cdf_mknam(c1,c2)
!  A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER*40 bsc_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      bsc_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      RETURN
       
      END FUNCTION bsc_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------

      END MODULE bsc_cdf
