      SUBROUTINE svdsolveb(m, n, mp, np, wcut, file, a, b, x)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, n, mp, np
      REAL(rprec) wcut
      REAL(rprec), DIMENSION(mp,np) :: a
      REAL(rprec), DIMENSION(n) :: b, x
      CHARACTER file*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n1, nonzero, nlast, istat, i, j, iunit
      REAL(rprec), ALLOCATABLE :: u(:,:), v(:,:), w(:)
      REAL(rprec) :: wmax, wcuta, udotb
C-----------------------------------------------

c       Solves Matrix equation A X = B for X using SVD method
c       Allows supression of weights, storage and reCALL of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu

c       Almost same as SvdSolveUV except this one stores basis vectors
c       rather than U, w, V. This one is a bit faster than SvdSolveUV, so
c       USE this one when you want to change ONLY weights, not B or A

c       Inputs:
c       A(M,N) - Matrix with physical dimensions (Mp,Np)
c       B(N)   - R.H.S. of A X = B
c       M      - Number of rows of A to use
c       N      - Number of columns of A to use
c       Mp,NP  - Storage dimensions of A(Mp,Np)
c       Wcut  - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wMAX < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = INTEGER > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c       Output:
c       X(N)   - Solution of A X = B
c
      n1 = n

c     ALLOCATE local arrays
      ALLOCATE (u(mp,np), w(np), v(np,np), stat=istat)

      IF (istat .eq. 0) THEN                       !have enough memory

         iunit = 41

         x(:n1) = zero                    !zero answer from previous CALL
         wcuta = ABS(wcut)
c...............................
c     If Wcut < 0 USE precalculated weights to get answer fast
c        USE wcuta = ABS(wcut) for everything
c     If ANY error happens WHILE READing "file", DO full svd calculation
         IF (wcut .lt. 0) THEN         !Read file from previous calculation
            CALL safe_open(iunit, istat, file, 'old', 'formatted')
            IF (istat .ne. 0) GOTO 98
            READ (iunit, *, err=98)
            READ (iunit, *, err=98)
            READ (iunit, *, err=98) m, n, nonzero
c                           !Read exactly the same way they were written
            DO i = 1, nonzero
               READ (iunit, *, err=98) w(i)         !READ i-th weight
               DO j = 1, n                 !READ i-th basis vector/ w(i)
c                             !READ j-th ENTRY (row) of i-th column of V
                  READ (iunit, *, err=98) v(j,i)
               END DO
            END DO
            CLOSE(unit=iunit)

c        Decide how many weights to keep
            IF (wcuta .gt. 1) nlast = wcuta         !Keep Nlast weights
            IF (wcuta .lt. 1) THEN                  !cutoff small weights
               wMAX = w(1)                          !weights are already ordered
               nlast = 0
               DO i = 1, nonzero
                  IF (w(i) .gt. wmax*wcuta) THEN
                     nlast = nlast + 1              !accept this weight
                  ELSE
                     GOTO 96
                  END IF
               END DO
            END IF

c        Just find linear SUM of Nlast vectors (weights were already included)
   96       CONTINUE
            IF (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
            DO i = 1, nlast             !add ith column of V(*,i) to X(*)
               x(:n) = x(:n) + v(:n,i)  !add to the j-th ENTRY of X
            END DO
            GOTO 999                   !ALL done
         ENDIF                          !End of precalculated branch Wcut < 0
   98    CONTINUE
c.......................................
c     If Wcut .ge. 0, DO full svd calculation
c     Initialize ALL to zero to wipe out effects of old CALL
         DO j = 1, n
            w(j) = zero                           !Zero ALL weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:m,j) = zero
         END DO

         CALL svdcmp (u, m, n, mp, np, w, v)     !Do SVD decomposition

c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         CALL sortsvd (m, n, mp, np, w, u, v)

c       Find the number of nonzero weights (already dcreasing ordered)
         DO nonzero = n, 1, -1
c                                   !Found first nonzero weight, get out
            IF (w(nonzero) .ne. 0) EXIT
         END DO
         IF (nonzero .gt. 0) THEN
            IF (wcuta .gt. 1) THEN

c        Decide how many weights to keep
               nlast = wcuta                     !Keep Nlast weights
            ELSE                                 !cutoff small weights
               wMAX = w(1)                  !weights are already ordered
               nlast = 0
               DO i = 1, nonzero
                  IF (w(i) .gt. wmax*wcuta) THEN
                     nlast = nlast + 1           !accept this weight
                  ELSE
                     GOTO 97
                  END IF
               END DO
            END IF

c       Find solution X and basis vectors (put into columns of V)
   97       CONTINUE
            IF (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
            DO i = 1, nlast  !Calc i-th coeff (U.b/w) for noNonzero w(i)
c        First calculate i-th coeff (U(i)*b/w(i)) in the SUM in eq 14.3.17 (NR)
c                                Dot product of i-th column of U with B
c                                summed over ALL plasma points
               udotb = SUM(u(:m,i)*b(:m))
c                        !=(Ui.b/wi) in eq 14.3.17, saves many divisions
               udotb = udotb/w(i)
c         Now run DOwn the i-th column of the vector V(j,i)
c                             j-th ENTRY (row) of vectors in eq 14.3.17
c                             Now V is the basis vector times weight
               v(:n,i) = udotb*v(:n,i)
c                             add it to the j-th ENTRY of X
               x(:n) = x(:n) + v(:n,i)
                                  !This did the svdbksb, faster than NR
            END DO                !and also calculated the basis vectors

c       Write all weights and weighted basis functions into SVD file (Wcut=0)
            IF (wcut .eq. 0) THEN
               CALL safe_open(iunit, istat, file, 'unknown','formatted')
               WRITE (iunit, *) 'Max w = ', w(1), ', Min w = ',
     1                          w(nonzero)
               WRITE (iunit, *) 'Ratio = ', w(1)/w(nonzero)
               WRITE (iunit, *) m, n, nonzero
               DO i = 1, nonzero             !WRITE in DECREASING order of weights
                  WRITE (iunit, *) w(i)      !WRITE i-th weight
                  DO j = 1, n                !WRITE i-th basis vector/ w(i)
                                             !WRITE j-th ENTRY (row) of i-th column of V
                     WRITE (iunit, *) v(j,i)
                  END DO
               END DO
               CLOSE(unit=iunit)
            END IF
c................................................
         END IF
      END IF
  999 CONTINUE

      DEALLOCATE (u, w, v, stat=istat)

      END SUBROUTINE svdsolveb
