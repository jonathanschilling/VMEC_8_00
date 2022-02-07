      SUBROUTINE svdinv(m, a, wcut, file)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m
      REAL(rprec) wcut
      CHARACTER file*(*)
      REAL(rprec), DIMENSION(m,m) :: a
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nonzero, nlast, istat, i, j, iunit
      REAL(rprec), ALLOCATABLE :: u(:,:), w(:), v(:,:)
      REAL(rprec) :: wmax, wcuta
C-----------------------------------------------

c       Inverts Square Matrix A using SVD method
c       Allows supression of weights, storage and reCALL of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu

c       A(M,M) - Matrix on input, Inverse on output
c       M      - Size of A
c       Wcut   - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wMAX < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = INTEGER > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c
c       Note: M is same as N, Mp, and Np in this routine

c     ALLOCATE local arrays
      ALLOCATE (u(m,m), w(m), v(m,m), stat=istat)
      IF (istat .eq. 0) THEN                       !have enough memory

         iunit = 41
         wcuta = ABS(wcut)
c...............................
c     If Wcut < 0 USE precalculated weights to get answer fast
c        USE wcuta = ABS(wcut) for everything
c     If ANY error happens while reading "file", DO full svd calculation
         IF (wcut .lt. 0) THEN         !Read file from previous calculation
            CALL safe_open(iunit, istat, file, 'old', 'formatted')
            IF (istat .ne. 0) GOTO 98
            READ (iunit, *, err=98)
            READ (iunit, *, err=98)
            READ (iunit, *, err=98) m, nonzero
c                           !Read exactly the same way they were written
            DO i = 1, nonzero
               READ (iunit, *, err=98) w(i)         !READ i-th weight
               DO j = 1, m
                  READ (iunit, *, err=98) v(j,i), u(j,i)
               END DO
            END DO
            CLOSE(unit=iunit)

c        Decide how many weights to keep
            nlast = wcuta                        !remember Wcut is < 0
            nlast = MIN(nonzero,nlast)      !Do not USE nonzero weights

            GOTO 99           !bypass svdcmp and GOTO matrix inversion
         ENDIF               !End of precalculated U,w,V branch Wcut < 0
   98    CONTINUE
c.......................................
c     If Wcut .ge. 0, DO full svd calculation and SAVE V,w,U
c     Initialize ALL to zero to wipe out effects of old CALL
         DO j = 1, m
            w(j) = zero                           !Zero ALL weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:m,j) = zero
         END DO

         CALL svdcmp (u, m, m, m, m, w, v)       !Do SVD decomposition

c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         CALL sortsvd (m, m, m, m, w, u, v)

c       Find the number of nonzero weights (already dcreasing ordered)
         DO nonzero = m, 1, -1
c                                   !Found first nonzero weight, get out
            IF (w(nonzero) .ne. 0) EXIT
         END DO
         IF (nonzero .le. 0) GOTO 999

c       Write ALL weights and U, V into SVD file (IF Wcut=0)
         IF (wcut .eq. 0) THEN    !no nonzero weights, you must be kidding
            CALL safe_open(iunit, istat, file, 'unknown', 'formatted')
            WRITE (iunit, *) 'Max w = ', w(1), ', Min w = ', w(nonzero)
            WRITE (iunit, *) 'Ratio = ', w(1)/w(nonzero)
            WRITE (iunit, *) m, nonzero
c                          !WRITE exactly the same way they were written
            DO i = 1, nonzero
               WRITE (iunit, *) w(i)                !WRITE i-th weight
               DO j = 1, m
                  WRITE (iunit, *) v(j,i), u(j,i)
               END DO
            END DO
            CLOSE(unit=iunit)
         END IF
c.......................................

c       Following gets done in both cases (using saved or recalculating)
c        Decide how many weights to keep
   99    CONTINUE
         IF (wcuta .gt. 1) THEN
            nlast = wcuta
         ELSE                                    !cutoff small weights
            wMAX = w(1)                     !weights are already ordered
            nlast = 0
            DO i = 1, nonzero
               IF (w(i) .gt. wmax*wcuta) THEN
                  nlast = nlast + 1              !accept this weight
               ELSE
                  GOTO 96
               END IF
            END DO
         END IF

c       Find Inverse of A = V [diag(1/wi)] Utranspose, RETURN in A
c       First DO [diag(1/wi)] Utranspose, store answer in U again
   96    CONTINUE
         IF (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
         DO i = 1, nlast
c                                         !divide ith row of Utr by w(i)
            u(:m,i) = u(:m,i)/w(i)
         END DO
c        Zero the infinite 1/weights
         DO i = nlast + 1, m
            u(:m,i) = zero                !multiply ith row of Utr by zero
         END DO
c       Next multiply by V matrix to get A inverse, put it in A
         DO i = 1, m                    !multiplt ith row vector of V by
            DO j = 1, m                         !jth column vector of wU
                                                !vector multiply V and wU
               a(i,j) = SUM(v(i,:m)*u(j,:m))    !Put inverse back in A
            END DO
         END DO
c................................................
      END IF
  999 CONTINUE

      DEALLOCATE (u, w, v, stat=istat)

      END SUBROUTINE svdinv
