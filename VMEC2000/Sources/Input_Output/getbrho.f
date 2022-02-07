      SUBROUTINE getbrho (brhomn, frho, bsupu, bsupv, mmax, nmax)
      USE stel_kinds
      USE vmec_input, ONLY: nfp, nzeta
      USE vmec_dim, ONLY: ntheta1, ntheta2, ntheta3
      USE vmec_persistent, ONLY: cosmu, sinmu, cosnv, sinnv
      USE vmec_main, ONLY: r0scale
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: mmax, nmax
      REAL(rprec), INTENT(out) :: brhomn(0:mmax, -nmax:nmax)
      REAL(rprec), DIMENSION(nzeta, ntheta3), INTENT(in) ::
     1    bsupu, bsupv, frho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, m, n, nmax1, itotal, ijtot, mntot
      REAL(rprec) :: ccmn, ssmn, dm, dn, termsc, termcs
      REAL(rprec), ALLOCATABLE :: amatrix(:,:), save_matrix(:,:),
     1   brhs(:)
      LOGICAL :: lpior0
      EXTERNAL solver
C-----------------------------------------------
!
!     Solves the radial force balance B dot Brho = Fs for Brho in real space using collocation
!     Here, Fs = frho(mn) is the Fourier transform SQRT(g)*F (part of radial force
!     balance sans the B dot Brho term)
!

      nmax1 = MAX(0,nmax-1)
      itotal = ntheta2*nzeta - 2*nmax1
      ALLOCATE (amatrix(itotal, itotal),
     1      brhs(itotal), save_matrix(itotal, itotal), stat=m)
      IF (m .ne. 0) STOP 'Allocation error in getbrho'


      amatrix = 0

!
!     BRHO = BSC(M,N)*SIN(MU)COS(NV) + BCS(M,N)*COS(MU)SIN(NV)
!          + BCC(M,N)*COS(MU)COS(NV) + BSS(M,N)*SIN(MU)SIN(NV)   (ONLY IF ASYMMETRIC MODE ALLOWED)
!

      ijtot = 0
      brhs = 0

      DO i = 1, ntheta2
         DO j = 1, nzeta
!           IGNORE u=0,pi POINTS FOR v > pi: REFLECTIONAL SYMMETRY
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            IF (lpior0) CYCLE
            ijtot = ijtot + 1
            brhs(ijtot) = frho(j,i)
            mntot = 0
            DO m = 0, mmax
               DO n = 0, nmax
                  ccmn = cosmu(i,m)*cosnv(j,n)
                  ssmn = sinmu(i,m)*sinnv(j,n)
                  dm = m * bsupu(j,i)
                  dn = n * bsupv(j,i) * nfp
                  termsc = dm*ccmn - dn*ssmn
                  termcs =-dm*ssmn + dn*ccmn
                  IF (n.eq.0 .or. n.eq.nmax) THEN
                     mntot = mntot + 1
                     IF (m .gt. 0) THEN
                        amatrix(ijtot,mntot) = termsc     !!ONLY bsc != 0 for n=0, nmax1
                     ELSE IF (n .eq. 0) THEN
                        amatrix(ijtot,mntot) = bsupv(j,i) !!pedestal for m=0,n=0 mode, which should = 0
                     ELSE
                        amatrix(ijtot,mntot) = termcs     !!bcs(m=0,n=nmax)
                     END IF
                  ELSE IF (m.eq.0 .or. m.eq.mmax) THEN
                     mntot = mntot + 1
                     amatrix(ijtot,mntot) = termcs        !!ONLY bcs != 0 for m=0,mmax
                  ELSE
                     amatrix(ijtot,mntot+1) = termsc
                     amatrix(ijtot,mntot+2) = termcs
                     mntot = mntot + 2
                  END IF
               END DO
            END DO
         END DO
      END DO

      save_matrix = amatrix

      IF (ijtot .ne. itotal .or. mntot .ne. itotal) THEN
         PRINT *,' itotal = ', itotal,' ijtot = ', ijtot,
     1   ' mntot = ', mntot
         STOP
      END IF
      IF (mmax+1 .ne. ntheta2) STOP 'Error 1 in getbrho'
      IF (nmax   .ne. nzeta/2) STOP 'Error 2 in getbrho'

      CALL solver (amatrix, brhs, itotal)

!
!     CHECK SOLUTION FROM SOLVER
!

!      ijtot = 0
!      DO i = 1, ntheta2
!         DO j = 1, nzeta
!           lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
!           IF (lpior0) CYCLE
!           ijtot = ijtot + 1
!           amn = SUM(save_matrix(ijtot,:)*brhs(:))
!           IF (ABS(frho(j,i) - amn) .gt. 1.e-8_dp*ABS(amn))
!    1      PRINT *,' i = ',i,' j = ',j,' Original force = ',
!    2      frho(j,i),' Final force = ', amn
!        END DO
!     END DO

!
!     CONVERT TO BS*SIN(MU - NV) REPRESENTATION
!
      mntot = 0
      brhomn = 0
      DO m = 0, mmax
         DO n = 0, nmax
            IF (n.eq.0 .or. n.eq.nmax) THEN
               mntot = mntot + 1
               IF (m .gt. 0) THEN
                  IF (n .eq. 0) THEN
                     brhomn(m,n) = brhs(mntot)         !!bcs(m,0), m > 0
                  ELSE
                     brhomn(m,n) = p5*brhs(mntot)      !!bsc(m,nmax), m > 0
                     brhomn(m,-n)= p5*brhs(mntot)
                  END IF
               ELSE IF (n .ne. 0) THEN
                  brhomn(m,n) = -brhs(mntot)            !!bcs(0,nmax): needed for derivative
               END IF
            ELSE IF (m.eq.0 .or. m.eq.mmax) THEN
               mntot = mntot + 1
               IF (m .eq. 0) THEN
                  brhomn(m,n)  = -brhs(mntot)          !!bcs(0,n) for n !=0, nmax
               ELSE
                  brhomn(m,n)  = -p5*brhs(mntot)       !!bcs(mmax,n)
                  brhomn(m,-n) =  p5*brhs(mntot)
               END IF
            ELSE
               brhomn(m,n)  = p5*(brhs(mntot+1) - brhs(mntot+2))
               brhomn(m,-n) = p5*(brhs(mntot+1) + brhs(mntot+2))
               mntot = mntot + 2
            END IF
         END DO
      END DO

      IF (mntot .ne. ijtot) STOP 'mntot != ijtot at END of getbrho'


      DEALLOCATE (amatrix, save_matrix, brhs)

      END SUBROUTINE getbrho
