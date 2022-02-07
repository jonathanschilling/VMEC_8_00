      PROGRAM vmec
      USE vmec_input
      USE vmec_seq
      USE safe_open_mod
      USE vparams, ONLY: nlog, nlog0
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nseq0 = 12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numargs, ierr_vmec, index_end,
     1   iopen, isnml, iread, iseq, index_seq,
     2   index_dat, ireset, iunit
      CHARACTER*120 :: input_file, seq_ext, reset_file_name, arg
      CHARACTER*120 :: log_file
      CHARACTER*120, DIMENSION(10) :: command_arg
      LOGICAL :: lfirst=.true., lreseta, lscreen
C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the PROGRAM VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report ANY problems or comments
!       to him.  As a BETA version, this PROGRAM is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       WHERE phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For ALL but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s FUNCTION table for the
!       vacuum magnetic field(s) and, IF data analysis is to be done, flux and
!       field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
!       called at an arbitrary spatial location and which should RETURN the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required IF
!       equilibrium reconstruction is to be done.)
!
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.
!
!
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to ANY one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these PARAMETER statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.
!
!
!       Added features since last edition
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!***

!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO

      lreseta = .true.            !!Default value: runvmec MUST be called this way the first time
      lscreen = .true.

      IF (numargs .lt. 1) THEN
         STOP 'Invalid command line'
      ELSE IF (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') THEN
         PRINT *,
     1   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         PRINT *
         PRINT *,' For example: '
         PRINT *,'    xvmec input.tftr OR xvmec tftr ',
     1           'OR xvmec ../input.tftr'
         PRINT *
         PRINT *,' Sequence files, containing a LIST of input files',
     1           ' are also allowed: '
         PRINT *,'    xvmec input.tftr_runs'
         PRINT *
         PRINT *,' Here, input.tftr_runs CONTAINS a &VSEQ NAMELIST',
     1           ' ENTRY'
         PRINT *
         PRINT *,' Additional (optional) command arguments are',
     1           ' allowed:'
         PRINT *
         PRINT *,'    xvmec <filename> noscreen F reset_wout_file'
         PRINT *
         PRINT *,' noscreen: supresses ALL output to screen ',
     1           ' (default, or "screen", displays output)'
         PRINT *,' F (or T): IF "T", forces reset on',
     1           ' a coarse mesh (used for sequencing control)'
         PRINT *,' name of reset wout file (defaults to this extension)'

         STOP
      ELSE IF (numargs .gt. 1) THEN
         arg = command_arg(2)
         IF (TRIM(arg).eq.'noscreen' .or. TRIM(arg).eq.'NOSCREEN')
     1      lscreen = .false.
      END IF
      IF (numargs .gt. 2) THEN
          arg = command_arg(3)
          IF (arg(1:1).eq.'f' .or. arg(1:1).eq.'F') lreseta = .false.
      END IF
      IF (numargs .gt. 3) THEN
          reset_file_name = command_arg(4)
      END IF


!
!     Determine type of file opened (sequential or input-data)
!     ARG1 (char var)
!          By DEFAULT, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to OPEN file ARG1 (full path + file name).
!                  Look for the VSEQ NAMELIST to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i]
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence NAMELIST, THEN the data file
!                  ARG1 or input.ARG1 is READ directly, with NSEQ=1.
!
      arg = command_arg(1)
      index_dat = index(arg,'.')
      index_end = len_trim(arg)
      IF (index_dat .gt. 0) THEN
         seq_ext  = arg(index_dat+1:index_end)
         input_file = TRIM(arg)
      ELSE
         seq_ext = TRIM(arg)
         input_file = 'input.'//TRIM(seq_ext)
      END IF

      IF (numargs .le. 3) reset_file_name = 'wout.' // seq_ext

      nseq = 1
      nseq_select(1) = 1
      extension(1) = input_file
!
!     READ IN NAMELIST VSEQ TO GET ARRAY
!     OF INPUT FILE EXTENSIONS AND INDEXING ARRAY, NSEQ_select
!
      nlog = nlog0
      iunit = nseq0
      DO iseq = 1, 2
         IF (iseq .eq. 1) THEN
           arg = input_file
         ELSE
           arg = seq_ext
         END IF
         CALL safe_open(iunit, iopen, TRIM(arg), 'old', 'formatted')
         IF (iopen .eq. 0) THEN
           CALL read_namelist (iunit, isnml, 'vseq')
           IF (isnml.eq.0 .and. nseq .gt. nseqmax) STOP 'NSEQ>NSEQMAX'

!
!       OPEN FILE FOR STORING SEQUENTIAL RUN HISTORY
!
           IF (isnml .eq. 0) THEN
              log_file = 'log.'//seq_ext

              CALL safe_open(nlog, iread, log_file, 'replace',
     1           'formatted')
              IF (iread .ne. 0) THEN
                 PRINT *, log_file,
     1           ' LOG FILE IS INACCESSIBLE: IOSTAT= ',iread
                 STOP 3
              ELSE
                 EXIT        !!Break out of loop
              END IF
           ENDIF
        ENDIF

        CLOSE (iunit)

      END DO

!
!     CALL EQUILIBRIUM SOLVER
!
!     nseq_select:      IF sequence file (VSEQ NAMELIST given with nseq >0)
!                       array giving indices into EXTENSION array prescribing
!                       the order in which the input files are run by VMEC
!     nseq:             number of sequential VMEC runs to make
!
!
!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
!
      DO iseq = 1, nseq
         index_seq = nseq_select(iseq)
         ireset = 0
         ierr_vmec = 0
         IF (iseq .gt. 1) reset_file_name =
     1       'wout.' // TRIM(extension(index_seq))
 100     CONTINUE
         CALL runvmec (extension(index_seq), iseq-1, lreseta, ierr_vmec,
     1                 ireset, lfirst, lscreen, reset_file_name)
         lfirst = .false.
         IF(ierr_vmec == 4)then
           IF(.not.lmoreiter) ierr_vmec = 0
         ENDIF
         SELECT CASE (ierr_vmec)
!        CASE (1:2)    !BAD JACOBIAN AFTER 75 ITERATIONS...
!          ireset = ireset + 1
!          lreseta = .true.
!          IF (ireset .le. 2) GOTO 100
         CASE (4)                                !Try a few more iterations
           ireset = ireset + 1
           lreseta = .false.
           IF (ireset .le. 1) THEN
              IF (lscreen) WRITE (6, '(/,1x,a)')
     1           'RUNNING A FEW MORE ITERATIONS THAN REQUESTED'
              GOTO 100
           ELSE IF (lscreen) THEN
              PRINT *, 'DECREASE DELT OR INCREASE NITER'
           ENDIF
         CASE (6)    !BAD JACOBIAN AFTER AXIS RESET: TRY DECREASING TO NS=3
           ireset = ireset + 1
           lreseta = .true.
           IF (ireset .le. 1) GOTO 100
         CASE DEFAULT
           lreseta = .false.
         END SELECT
      END DO

!
!     FREE ANY LONG-TERM (PERSISTENT THROUGH ISEQ > 1, OR XC, SCALXC FOR
!     ITERATIVE OPTIMIZATION) POINTERS
!
      CALL free_persistent_mem

      CLOSE (nlog)

      END PROGRAM vmec
