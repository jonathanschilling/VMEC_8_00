      SUBROUTINE totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1,
     1   lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
C-----------------------------------------------
      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncs,2)
      joff = LBOUND(rmncs,3)

      z00b = zmncc(ns,ioff,joff)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC totzspa'

!
!                INITIALIZATION BLOCK
!
      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0

      IF (jlam(m0) .gt. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work1(j1l:nsl,1) = work1(j1l:nsl,1)
     1                          + rmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,6) = work1(j1l:nsl,6) 
     1                          + zmncc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,10) = work1(j1l:nsl,10)
     1                          + lmncc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               work1(j1l:nsl,2) = work1(j1l:nsl,2)
     1                          + rmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7)
     1                          + zmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12) 
     1                          + lmncc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO

!
!        INVERSE TRANSFORM IN M-THETA
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux

            IF (.not.lthreed) CYCLE
               
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               sinmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               cosmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               cosmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               sinmu(i,m) + work1(:,4)*cosmu(i,m)
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               sinmu(i,m) + work1(:,8)*cosmu(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               sinmu(i,m)+work1(:,12)*cosmu(i,m))
         END DO
      END DO

      DEALLOCATE (work1)

      END SUBROUTINE totzspa
