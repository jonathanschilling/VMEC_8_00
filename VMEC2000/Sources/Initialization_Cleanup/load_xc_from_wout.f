      SUBROUTINE load_xc_from_wout(rmn, zmn, lmn, ntor_in, mpol1_in,
     1    ns_in, reset_file, lreset)
      USE read_wout_mod, ONLY: rmnc, zmns, lmns, rmns, zmnc, lmnc,
     1    xm, xn, ntor, ns,
     2    nfp, mnmax, read_wout_file, READ_wout_deallocate
      USE vmec_params, ONLY: rprec => rprec1, mscale, nscale, ntmax
      USE vmec_dim, ONLY: mpol1
      USE vparams, ONLY: one, zero
      USE vmec_input, ONLY: lasym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ns_in, mpol1_in, ntor_in
      REAL(rprec), DIMENSION(ns_in,0:ntor_in,0:mpol1_in,ntmax),
     1   INTENT(out) :: rmn, zmn, lmn
      CHARACTER*(*) :: reset_file
      LOGICAL lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: rbcc=1, rbss=2, rbcs=3, rbsc=4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ierr, mn, m, n, n1, zbcs, zbsc, zbcc, zbss
      REAL(rprec) :: t1, t2
C-----------------------------------------------


!
!     THIS ALLOWS SEQUENTIAL RUNNING OF VMEC FROM THE COMMAND LINE
!     i.e., WHEN VMEC INTERNAL ARRAYS ARE NOT KEPT IN MEMORY (compared to sequence file input)
!     THIS IS THE CASE WHEN VMEC IS CALLED FROM, SAY, THE OPTIMIZATION CODE
!
      CALL read_wout_file (reset_file(5:), ierr)

      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening/reading wout file in VMEC load_xc!'
         lreset = .true.
         RETURN
      END IF

      lreset = .false.

      IF (ns_in .ne. ns) THEN
         lreset = .true.
         PRINT *, 'ns_in != ns in load_xc'
         RETURN
      END IF

      IF (ntor_in  .ne. ntor ) STOP 'ntor_in != ntor in load_xc'
      IF (mpol1_in .ne. mpol1) STOP 'mpol1_in != mpol1 in load_xc'
      IF (nfp .eq. 0) STOP 'nfp = 0 in load_xc'

      rmn = zero
      zmn = zero
      lmn = zero

      zbcs = rbcc
      zbsc = rbss
      zbcc = rbcs
      zbss = rbsc

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))/nfp
         n1 = ABS(n)
         t1 = one/(mscale(m)*nscale(n1))
         t2 = t1
         IF (n .lt. 0) t2 = -t2
         IF (n .eq. 0) t2 = zero
         rmn(:ns, n1, m, rbcc) = rmn(:ns, n1, m, rbcc) + t1*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbsc) = zmn(:ns, n1, m, zbsc) + t1*zmns(mn,:ns)
         lmn(:ns, n1, m, zbsc) = lmn(:ns, n1, m, zbsc) + t1*lmns(mn,:ns)
         rmn(:ns, n1, m, rbss) = rmn(:ns, n1, m, rbss) + t2*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbcs) = zmn(:ns, n1, m, zbcs) - t2*zmns(mn,:ns)
         lmn(:ns, n1, m, zbcs) = lmn(:ns, n1, m, zbcs) - t2*lmns(mn,:ns)
         IF (lasym) THEN
         rmn(:ns, n1, m, rbcs) = rmn(:ns, n1, m, rbcs) - t2*rmns(mn,:ns)
         zmn(:ns, n1, m, zbcc) = zmn(:ns, n1, m, zbcc) + t1*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbcc) = lmn(:ns, n1, m, zbcc) + t1*lmnc(mn,:ns)
         rmn(:ns, n1, m, rbsc) = rmn(:ns, n1, m, rbsc) + t1*rmns(mn,:ns)
         zmn(:ns, n1, m, zbss) = zmn(:ns, n1, m, zbss) + t2*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbss) = lmn(:ns, n1, m, zbss) + t2*lmnc(mn,:ns)
         END IF
         IF (m .eq. 0) THEN
            rmn(:ns, n1, m, rbss) = zero
            zmn(:ns, n1, m, zbsc) = zero
            lmn(:ns, n1, m, zbsc) = zero
            IF (lasym) THEN
            rmn(:ns, n1, m, rbsc) = zero
            zmn(:ns, n1, m, zbss) = zero
            lmn(:ns, n1, m, zbss) = zero
            END IF
         END IF
      END DO


      CALL read_wout_deallocate

      END SUBROUTINE load_xc_from_wout
