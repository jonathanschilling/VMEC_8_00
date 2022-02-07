      SUBROUTINE fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag,
     1    ncnt, epsfcn, wa, time, fnorm_min, x_min, fvec_min)
      USE stel_kinds
!DEC$ IF .NOT.DEFINED (MPI_OPT)
      USE fdjac_mod, m1=>m, n1=>n, eps1=>eps, ncnt1=>ncnt
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, ldfjac, ncnt
      INTEGER, TARGET :: iflag
      REAL(rprec) ::  epsfcn, time
      REAL(rprec), DIMENSION(n), TARGET :: x(n), wa(m)
      REAL(rprec), DIMENSION(m), INTENT(in) :: fvec
      REAL(rprec), DIMENSION(ldfjac,n), INTENT(out) :: fjac
      REAL(rprec), INTENT(out) :: fnorm_min, x_min(n), fvec_min(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max
      REAL(rprec) :: eps, epsmch, h, dpmpar, temp, cur_norm
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn, dpmpar, multiprocess, fdjac_parallel
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     Here
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an EXTernal statement in the user calling
c         program (see LMDIF1 for documentation), and should be written as follows:
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag,ncnt)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c       fjac is an output m by n array which CONTAINS the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. IF epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... ABS,max,sqrt
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))
!
!     Load MODULE values. Pointers will automatically update TARGETs...
!     Prepare for multi-processing...
!
      m1 = m
      n1 = n
      ncnt1 = ncnt
      eps1 = eps
      xp => x
      wap => wa

!     Find MIN chisq = fnorm**2 state for this jacobian evaluation
!     (Do NOT retain from previous evaluations, or could get into a non-converging loop...)

      fnorm_min = HUGE(fnorm_min)

      CALL system_clock(ic1, irate)
      CALL multiprocess(n, max_processors, fdjac_parallel, fcn)
      CALL system_clock(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

c     DO j = 1, n
c         temp = x(j)
c         h = eps*ABS(temp)
c         IF (h .eq. zero) h = eps
c         x(j) = temp + h
c         CALL fcn (m, n, x, wa, iflag, ncnt)
c         IF (iflag .lt. 0) EXIT
c         x(j) = temp
c         fjac(:m,j) = (wa - fvec)/h
c      END DO


      cur_norm = enorm(m,fvec)      ! WHERE are we now?

      DO j = 1, n

        READ (j+1000, iostat=iread) istat, iflag, h, temp
        IF (iREAD .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1     ' in fdjac2: IOSTAT = ', iread
           iflag = -14
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *) 'Wrong value for INDEX j READ in fdjac2'
           iflag = -14
        END IF

        IF (iflag .ne. 0) EXIT

!DEC$ IF DEFINED (CRAY)
        DO k = 1, m
           READ (j+1000) wa(k)
        END DO
!DEC$ ELSE
        READ (j+1000) wa
!DEC$ ENDIF
        fjac(:m,j) = (wa - fvec)/h

        IF( temp > cur_norm) flip(j) = .not. flip(j)  ! flip for next time

        IF( temp < fnorm_min) THEN
           fnorm_min = temp
           fvec_min = wa
!DEC$ IF DEFINED (CRAY)
           DO k = 1, n
              READ (j+1000) x_min(k)
           END DO
!DEC$ ELSE
           READ (j+1000) x_min
!DEC$ ENDIF
        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

!
!     Do ANY special cleanup now for IFLAG = flag_cleanup
!
      iflag = flag_cleanup
      CALL fcn(m, n, x, wa, iflag, ncnt)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
!DEC$ ENDIF
      END SUBROUTINE fdjac2
