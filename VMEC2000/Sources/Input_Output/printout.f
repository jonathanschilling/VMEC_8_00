      SUBROUTINE printout(i0, delt0, w0, lscreen)
      USE vmec_main
      USE realspace
      USE vsvd
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: i0
      REAL(rprec) :: delt0, w0
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER*(*), PARAMETER :: 
     1   iter_line = " ITER    FSQR      FSQZ      FSQL   ",
     2   fsq_line  = "   fsqr      fsqz      fsql      DELT    ",
     3   raxis_line = " RAX(v=0) ",
     4   zaxis_line = "  ZAX(v=0)      "
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: betav, w, avm, den
      CHARACTER(len=LEN(iter_line) + LEN(fsq_line) +
     1          LEN(raxis_line) + LEN(zaxis_line)) :: 
     2          print_line
C-----------------------------------------------
      betav = wp/wb
      w = w0*twopi*twopi
      den = zero
      specw(1) = one
      CALL spectrum (xc(:irzloff), xc(1+irzloff:2*irzloff))
      den = SUM(phip(2:ns))
      avm = DOT_PRODUCT(phip(2:ns),specw(2:ns)+specw(:ns-1))
      avm = 0.5_dp*avm/den
      IF (ivac .ge. 1 .and. iter2.gt.1) delbsq =
     1     SUM(dbsq(:nznt)*wint(2:nrzt:ns))/
     1     SUM(bsqsav(:nznt,3)*wint(2:nrzt:ns))
      IF (i0.eq.1 .and. lfreeb) THEN
         print_line = iter_line // " " // raxis_line 
         IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
         IF (lscreen) PRINT 20, TRIM(print_line)
         print_line = iter_line // fsq_line // raxis_line
         IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
         IF (imatch_phiedge .eq. 1) THEN
            WRITE (nthreed, 15) TRIM(print_line)
         ELSE
            WRITE (nthreed, 16) TRIM(print_line)
         ENDIF
      ELSE IF (i0.eq.1 .and. .not.lfreeb) THEN
         print_line = raxis_line 
         IF (lasym) print_line = raxis_line // zaxis_line
         IF (lscreen) PRINT 30, iter_line, TRIM(print_line)
         print_line = iter_line // fsq_line // raxis_line // "     "
         IF (lasym) print_line = iter_line // fsq_line // raxis_line
     1                        // zaxis_line
         WRITE (nthreed, 25) TRIM(print_line)
      ENDIF
   15 FORMAT(/,a,6x,'WMHD      BETA      <M>   DEL-BSQ   FEDGE',/)
   16 FORMAT(/,a,6x,'WMHD      BETA     PHIEDGE  DEL-BSQ    FEDGE',/)
   20 FORMAT(/,a,6x,'WMHD      DEL-BSQ',/)
   25 FORMAT(/,a,6x,'WMHD      BETA      <M>        ',/)
   30 FORMAT(/,a,1x,a,5x,'WMHD',/)
      IF (.not. lasym) THEN
         IF (.not.lfreeb) THEN
            IF (lscreen) PRINT 45, i0, fsqr, fsqz, fsql, r00, w
            WRITE (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, w, betav, avm
         RETURN
         ENDIF
         IF (lscreen) PRINT 50, i0, fsqr, fsqz, fsql, r00, w, delbsq
         IF (imatch_phiedge .eq. 1) THEN
            WRITE (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, w, betav, avm, delbsq, fedge
         ELSE
            WRITE (nthreed, 42) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, w, betav, ABS(phiedge*phifac), delbsq, 
     2      fedge
         ENDIF
      
      ELSE
      
         IF (.not.lfreeb) THEN
            IF (lscreen) PRINT 65, i0, fsqr, fsqz, fsql, r00, z00, w
            WRITE (nthreed, 60) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, z00, w, betav, avm
         RETURN
         ENDIF
         IF (lscreen) PRINT 70, i0, fsqr, fsqz, fsql, r00, z00,w,delbsq
         IF (imatch_phiedge .eq. 1) THEN
            WRITE (nthreed, 60) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, z00, w, betav, avm, delbsq, fedge
         ELSE
            WRITE (nthreed, 60) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, 
     1      fsql1, delt0, r00, z00, w, betav, ABS(phiedge*phifac),
     2      delbsq,fedge
         ENDIF
      END IF

   40 FORMAT(i5,1p,7e10.2,e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
   42 FORMAT(i5,1p,7e10.2,e11.3,e12.4,2e11.3,0p,f7.3,1p,e9.2)
   45 FORMAT(i5,1p,3e10.2,e11.3,e12.4)
   50 FORMAT(i5,1p,3e10.2,e11.3,e12.4,e11.3)
   60 FORMAT(i5,1p,7e10.2,2e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
   65 FORMAT(i5,1p,3e10.2,2e11.3,e12.4)
   70 FORMAT(i5,1p,3e10.2,2e11.3,e12.4,e11.3)

!     IF (lrecon) THEN
!        IF (ystark(1)*ystark(isnodes).lt.zero .or.
!    1       ystark(2)*ystark(isnodes-1).lt.zero) THEN
!           IF (lscreen) PRINT 200
!           WRITE (nthreed, 200)
!        ENDIF
! 200    FORMAT(' CHECK THE SIGN OF MSE DATA: INCONSISTENT WITH ',
!    1      'CURRENT, TOROIDAL FIELD!')
!       IF (pressum0 .ne. zero .and. apres .ne. zero .and. lscreen)
!    1     PRINT 60, errsvd, pfac, phifac,
!    2               ABS(fsqsum0/pressum0), aminor/apres
!  60    FORMAT(2x,'chisq-err   = ',f8.3,3x,'pres-scale = ',f8.3,3x,
!    1      'phifac    = ',f8.3,/,2x,'<F00> force = ',f8.3,3x,
!    2      '<a>/apres  = ',f8.3)
!        IF (imatch_phiedge.eq.2 .and. ivac>1) THEN
!           IF (lscreen) PRINT 65, dlim_min, rlim_min, zlim_min
!  65       FORMAT(2x,'dlim (min)  = ',f8.3,3x,'at rlim    = ',f8.3,
!    1         '   zlim = ',f8.3/)
!        ELSE
!           IF (lscreen) PRINT *, ' '
!        ENDIF
!     ENDIF

      END SUBROUTINE printout
