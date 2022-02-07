      SUBROUTINE fileout(iseq, ier_flag, ireset, lscreen)
      USE vmec_main
      USE vac_persistent
      USE realspace
      USE vmec_params, ONLY: mscale, nscale, signgs, uminus,
     1      norm_term_flag, more_iter_flag, bad_jacobian_flag
      USE vforces
      USE vsvd
      USE xstuff, ONLY: xc, gc, xcdot
      USE precon2d, ONLY: lprec2d, rzl_save
      USE timer_sub

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iseq, ier_flag, ireset
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      LOGICAL, PARAMETER :: lreset_xc = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, istat1=0
      REAL(rprec), DIMENSION(:), POINTER :: lu, lv
      REAL(rprec), ALLOCATABLE :: br(:), bz(:)
      CHARACTER*(*), PARAMETER, DIMENSION(0:10) :: werror = (/
     1   'EXECUTION TERMINATED NORMALLY                            ',
     1   'INITIAL JACOBIAN CHANGED SIGN (IMPROVE INITIAL GUESS)    ',
     2   'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)         ',
     3   'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.            ',
     4   'FORCE RESIDUALS EXCEED FTOL: INCREASING NUMBER ITERATIONS',
     5   'ERROR READING INPUT FILE OR NAMELIST                     ',
     6   'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN        ',
     7   'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE              ',
     8   'NS ARRAY MUST NOT BE ALL ZEROES                          ',
     9   'ERROR READING MGRID FILE                                 ',
     A   'VAC-VMEC I_TF MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL  ' /)
      CHARACTER*(*), PARAMETER ::
     1    Warning = " Global memory deallocation in FILEOUT, error "
      LOGICAL :: log_open, lwrite
C-----------------------------------------------
   
      lu => czmn;   lv => crmn

!
!     COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
!     CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
!     AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
!
      iequi = 1
      lwrite = (ier_flag.eq.norm_term_flag) .or. 
     1         (ier_flag.eq.more_iter_flag .and. ireset.ge.1)

      IF (lwrite) THEN
!
!     The sign of the jacobian MUST multiply phi to get the physically
!     correct toroidal flux
!
         phi(1) = zero
         DO js = 2, ns
            phi(js) = phi(js-1) + phip(js)
         END DO
         phi = (signgs*twopi*hs)*phi

         CALL funct3d (lscreen, ier_flag)

         CALL second0 (teqfon)
         ALLOCATE(br(nrzt), bz(nrzt))
         CALL eqfor (br, bz, clmn, blmn, 
     1               rcon(1,1), xc)
         CALL second0 (teqfoff)
         timer(teqf) = timer(teqf) + teqfoff - teqfon
      END IF
  
!
!     MUST call WROUT to write error, if nothing else
!
      CALL second0 (twouton)

!DEC$ IF DEFINED (NETCDF)
      IF (laddout) CALL wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o,
     1                   czmn_e, crmn_e, azmn_e, ier_flag, lwrite)
      CALL wrout_cdf (bzmn_o, azmn_o, clmn, blmn, crmn_o, czmn_e,
     1     crmn_e, xcdot, gc, ier_flag, lwrite)
!DEC$ ELSE
      CALL wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o,
     1      czmn_e, crmn_e, azmn_e, ier_flag, lwrite)
!DEC$ ENDIF

      CALL second0 (twoutoff)
     
!
!     TESTING READ_WOUT MODULE WRITING ROUTINES
!
!      IF (lwrite) THEN
!         CALL TestWout(xc, br, bz, crmn_e, czmn_e)
!         DEALLOCATE (br, bz)
!      END IF

!     END TEST

      timer(twout) = timer(twout) + twoutoff - twouton
      timer(tsum)  = timer(tsum) + timer(twout) + timer(teqf)

      IF (lscreen) PRINT 10, ijacob
      IF (lscreen) PRINT 120, TRIM(werror(ier_flag)), input_extension
      IF (ier_flag .ne. norm_term_flag) GOTO 1000

      IF (nthreed .gt. 0) THEN
         WRITE (nthreed, 10) ijacob
         CALL write_times(nthreed, lscreen, lfreeb, lrecon, lprec2d)
         WRITE (nthreed, 120) TRIM(werror(ier_flag)), input_extension
      END IF
   10 FORMAT(/,'  NUMBER OF JACOBIAN RESETS = ',i4,/)
  120 FORMAT(/,2x,a,/,2x,'FILE : ',a/)
!
!     WRITE SEQUENCE HISTORY FILE
!
      INQUIRE(unit=nlog,opened=log_open)
      IF (lrecon .and. log_open) THEN
         IF (iseq .eq. 0) WRITE (nlog, 100)
         WRITE (nlog, 110) iseq + 1, iter2, total_chi_square_n,
     1      1.e-6_dp*ctor/mu0, 1.e-3_dp*ppeak/mu0, torflux, r00,
     2      timer(tsum), input_extension
      ENDIF

  100 FORMAT(' SEQ ITERS  CHISQ/N',
     1   '  TORCUR  PRESMAX  PHIEDGE     R00 CPU-TIME  EXTENSION')
  110 FORMAT(i4,i6,f8.2,3f9.2,f8.2,f9.2,2x,a20)

 1000 CONTINUE

!
!     DEALLOCATE GLOBAL MEMORY
!
      IF (ALLOCATED(cosmu))
     1  DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
     2  sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn, 
     3  cos01, sin01, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#1"

      IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, ixm, xm_nyq, xn_nyq, 
     1   jmin3, mscale, nscale, uminus, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#2"

      IF (ALLOCATED(tanu))
     1  DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv,
     2  sinu, cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
     3  cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#3"

      CALL free_mem_funct3d
      CALL free_mem_ns (lreset_xc)
      CALL free_mem_nunv
      IF (ALLOCATED(rzl_save)) DEALLOCATE(rzl_save)

!
!     CLOSE OPENED FILES
!
      IF (ier_flag .ne. more_iter_flag) CALL close_all_files

      END SUBROUTINE fileout
