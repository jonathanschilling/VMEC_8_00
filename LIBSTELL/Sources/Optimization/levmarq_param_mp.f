      SUBROUTINE levmarq_param_mp(x, wa1, wa2, wa3, wa4,
     1     nfev, m, n, iflag, fcn)
      USE lmpar_mod
      USE fdjac_mod, ONLY: flag_cleanup
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: n, m
      INTEGER :: iflag, nfev
      REAL(rprec), INTENT(in) :: x(n)
      REAL(rprec) :: wa1(n), wa2(n), wa3(n), wa4(m)
      EXTERNAL fcn
!DEC$ IF DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec), DIMENSION(11), PARAMETER :: factors =
     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 0.75_dp,
     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp, 2.1_dp /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iproc, iproc_min, nfact, num_lev, istat, ierr, j
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iflag_array
      REAL(rprec) :: scale_factor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fnorm_array
      CHARACTER*1 :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
!     Perform Numprocs different function calls (in parallel) find the minimum norm (chi-sq).
!     MPI calls are used to determine which processor has the minimum norm and then the
!     information is sent to all processors using MPI Broadcasts (modifications made by D. A. Spong 8/27/2000).

      nfact = SIZE(factors)
      num_lev = MIN(numprocs, n)
      ALLOCATE (iflag_array(numprocs), fnorm_array(numprocs),
     1          stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in levmarq_param_mp'

      iproc = myid + 1
      IF (lfirst_lm .and. num_lev > 2) THEN
!
!       Do an exponential spread the first time to see where we are
!
!SPH     scale_factor = EXP((iproc-1)*log(spread_ratio)/num_lev)
         scale_factor = 10._dp**(1._dp - iproc)
      ELSE IF (num_lev > 2*nfact) THEN
         scale_factor = (iproc*MAXVAL(factors))/num_lev
      ELSE IF (iproc .le. nfact) THEN
         scale_factor = factors(iproc)
      ELSE
         scale_factor = ((iproc-nfact)*MINVAL(factors))/
     1                   (num_lev-nfact)
      END IF

      delta = scale_factor*delta

      CALL lmpar (n, fjac, ldfjac, ipvt, diag, qtf,
     1            delta, par, wa1, wa2, wa3, wa4)
!
!     store the direction p and x + p. calculate the norm of p.
!
      IF (par.eq.zero .and. myid.ne.master) wa1 = wa1*scale_factor
      wa1 = -wa1
      wa2 = x + wa1
      wa3 = diag*wa1
      pnorm = enorm(n,wa3)
!
!     evaluate the function at x + p and calculate its norm.
!     Only DO for 0 <= myid < n processors to avoid clean-up problems (in lsfun1)
!
      IF (iproc .le. num_lev) THEN
         iflag = iproc
         CALL fcn (m, n, wa2, wa4, iflag, nfev)
      ELSE
         iflag = 0
      END IF

!
!     Gather iflag information to all processors and check for iflag < 0
!
      CALL MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_WORLD, ierr)
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in LMDIF'

      iflag = MINVAL(iflag_array)
      IF (iflag .lt. 0) RETURN

      IF (iproc .le. num_lev) THEN
         fnorm1 = enorm(m,wa4)
      ELSE
         fnorm1 = HUGE(fnorm1)
      END IF

!
!     Find processor with minimum fnorm1 value
!
      CALL MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in LMDIF'
      iflag_array = MINLOC(fnorm_array)
      iproc_min = iflag_array(1) - 1

      ext = ' '
      low_mark = ' '
      IF (myid .eq. master) ext = '*'
      IF (iproc .eq. iproc_min+1) low_mark = '*'
      IF (iproc .le. num_lev)
     1   WRITE(6, '(2x,i6,8x,i3,4x,2(3x,1es12.4,a),3x,1es12.4)')
     2         iproc+nfev, iproc, fnorm1**2, low_mark, par, ext,
     3         delta

      fnorm1 = MINVAL(fnorm_array)

      DEALLOCATE(iflag_array, fnorm_array)

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     ALL processors. wa3 is overwritten...
!
      CALL MPI_BCAST(pnorm,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(par,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(delta,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa1,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      IF (ierr .ne. 0) GOTO 3000

!
!     BROADCAST JACOBIAN fjac(:n,j), j=1,n ROW BY ROW (CHANGED IN LMPAR) TO OTHER PROCESSES
!
      DO j = 1, n
         IF (myid .eq. iproc_min) wa3(:n) = fjac(:n,j)
         CALL MPI_BCAST(wa3, n, MPI_REAL8, iproc_min,
     1        MPI_COMM_WORLD, ierr)
         IF (ierr .ne. 0) GOTO 3000
         IF (myid .ne. iproc_min) fjac(:n,j) = wa3(:n)
      END DO

!
!     CLEANUP AFTER LEVENBERG-MARQUARDT LOOP AS NEEDED (WA4 IS NOT CHANGED)
!
      iflag = flag_cleanup
      CALL fcn (m, n, x, wa4, iflag, nfev)             !Contains Bcast Barrier

      nfev = nfev + num_lev

      RETURN

 3000 CONTINUE

      WRITE (6, *) 'MPI_BCAST error in LEVMARQ_PARAM_MP, ierr = ', ierr
!DEC$ ENDIF
      END SUBROUTINE levmarq_param_mp
