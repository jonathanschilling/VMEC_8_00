      SUBROUTINE tomnspa(frzl_array, armn, brmn, crmn, azmn, bzmn,
     1   czmn, blmn, clmn, arcon, azcon)
      USE vmec_main
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: frzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) ::
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           frcs, frsc, fzcc, fzss, flcc, flss
      REAL(rprec), DIMENSION(ns*nzeta) :: temp1, temp3
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
C-----------------------------------------------
      frsc => frzl_array(:,:,:,rsc)               !!R-SIN(mu) COS(nv)
      fzcc => frzl_array(:,:,:,zcc+ntmax)         !!Z-COS(mu) COS(nv)
      flcc => frzl_array(:,:,:,zcc+2*ntmax)       !!L-COS(mu) COS(nv)
      IF (lthreed) THEN
         frcs => frzl_array(:,:,:,rcs)            !!R-COS(mu) SIN(nv)
         fzss => frzl_array(:,:,:,zss+ntmax)      !!Z-SIN(mu) SIN(nv)
         flss => frzl_array(:,:,:,zss+2*ntmax)    !!L-SIN(mu) SIN(nv)
      END IF

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC tomnspa'

      ioff = LBOUND(frcs,2)
      joff = LBOUND(frcs,3)

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1
!
!     BEGIN FOURIER TRANSFORM
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) TRANSFORM FIRST
!
         DO i = 1, ntheta2
            temp1(:) = armn(:,i,mparity) + xmpq(m,1)*arcon(:,i,mparity)
            temp3(:) = azmn(:,i,mparity) + xmpq(m,1)*azcon(:,i,mparity)
            work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m) 
     1                              + brmn(:,i,mparity)*cosmumi(i,m)
            work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m)
     1                              + bzmn(:,i,mparity)*sinmumi(i,m)
            work1(:,9) = work1(:,9) + blmn(:,i,mparity)*sinmumi(i,m)

            IF (.not.lthreed) CYCLE

            work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m)
     1                              + brmn(:,i,mparity)*sinmumi(i,m)
            work1(:,2) = work1(:,2) - crmn(:,i,mparity)*cosmui(i,m)
            work1(:,4) = work1(:,4) - crmn(:,i,mparity)*sinmui(i,m)
            work1(:,6) = work1(:,6) - czmn(:,i,mparity)*cosmui(i,m)
            work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m)
     1                              + bzmn(:,i,mparity)*cosmumi(i,m)
            work1(:,8) = work1(:,8) - czmn(:,i,mparity)*sinmui(i,m)
            work1(:,10) = work1(:,10) - clmn(:,i,mparity)*cosmui(i,m)
            work1(:,11) = work1(:,11) + blmn(:,i,mparity)*cosmumi(i,m)
            work1(:,12) = work1(:,12) - clmn(:,i,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,3)*cosnv(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,5)*cosnv(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,9)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,4)*sinnvn(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,6)*sinnvn(k,n)
               frcs(j2:jmax,ni,mj) = frcs(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,1)*sinnv(k,n) 
     2                             + work1(j2l:jmaxl,2)*cosnvn(k,n)
               fzss(j2:jmax,ni,mj) = fzss(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,7)*sinnv(k,n)
     2                             + work1(j2l:jmaxl,8)*cosnvn(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,10)*sinnvn(k,n)
               flss(jl:ns,ni,mj) = flss(jl:ns,ni,mj)
     1                           + work1(jll:nsl,11)*sinnv(k,n)
     2                           + work1(jll:nsl,12)*cosnvn(k,n)
            END DO
         END DO
      END DO

      DEALLOCATE (work1)

      END SUBROUTINE tomnspa
