      SUBROUTINE levmarq_param(x, wa1, wa2, wa3, wa4,
     1     time, nfev, m, n, iflag, fcn)
      USE fdjac_mod, mp => m, np => n
      USE lmpar_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nfev, n, m, iflag
      REAL(rprec) :: time
      REAL(rprec), TARGET :: x(n), wa1(n), wa2(n), wa3(n),
     1    wa4(m)
      EXTERNAL fcn
!DEC$ IF .NOT.DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max,
     1     jmin
      REAL(rprec), DIMENSION(num_lm_params) ::
     1      fnorm_min, pnorm_min, delta_min, par_min
      CHARACTER*1 :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL multiprocess, lmpar_parallel
C-----------------------------------------------
!
!     Initialize MODULE variables for USE in lmpar_parallel
!

      xp => x
      wap => wa1
      wa2p => wa2
      wa3p => wa3
      wa4p => wa4
      np = n
      mp = m
      ncnt = nfev

      CALL system_clock(ic1, irate)

      CALL multiprocess(num_lm_params, max_processors,
     1    lmpar_parallel, fcn)

      CALL system_clock(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

      nfev = nfev + num_lm_params

!
!     Read in optimal wa1, wa2, wa4, par, delta, pnorm, fnorm1 value from file
!
      DO j = 1, num_lm_params

        READ (j+1000, iostat=iread) istat, iflag, pnorm_min(j),
     1        fnorm_min(j), par_min(j), delta_min(j)
        IF (iREAD .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1       ' in levmarq_param', ' IOSTAT = ', iread
           iflag = -15
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *)
     1        'Incorrect value READ in for INDEX j in levmarq_param'
           iflag = -15
        END IF

        IF (iflag .ne. 0) RETURN

        IF (j .eq. 1) fnorm1 = fnorm_min(j)
        IF (fnorm_min(j) .le. fnorm1) THEN
           jmin = j
           fnorm1 = fnorm_min(jmin)
           pnorm  = pnorm_min(jmin)
           par    = par_min(jmin)
           delta  = delta_min(jmin)
!DEC$ IF DEFINED (CRAY)
           DO k = 1, n
              READ (j+1000) wa1(k), wa2(k)
              DO istat = 1, n
                 READ (j+1000) fjac(k, istat)
              END DO
           END DO
           DO k = 1, m
              READ (j+1000) wa4(k)
           END DO
!DEC$ ELSE
           READ (j+1000) wa1, wa2, wa4, fjac(1:n, 1:n)
!DEC$ ENDIF

        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

      DO j = 1, num_lm_params
         ext = ' '
         low_mark = ' '
         IF (j .eq. 1) ext = '*'
         IF (j .eq. jmin) low_mark = '*'
         WRITE (6, '(2x,i6,4x,2(3x,1es12.4,a),3x,1es12.4)') j+ncnt,
     1         fnorm_min(j)**2, low_mark, par_min(j), ext, delta_min(j)
      END DO

!
!     Do ANY special cleanup now for IFLAG = flag_cleanup. WA4 LEFT UNCHANGED
!
      iflag = flag_cleanup
      CALL fcn(m, n, x, wa4, iflag, ncnt)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
!DEC$ ENDIF
      END SUBROUTINE levmarq_param
