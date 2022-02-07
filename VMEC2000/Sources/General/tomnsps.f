      SUBROUTINE tomnsps(frzl_array, armn, brmn, crmn, azmn, 
     1   bzmn, czmn, blmn, clmn, arcon, azcon)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs,
     1                       mscale, nscale
      USE fbal, ONLY: lforbal, rru_fac, rzu_fac, frcc_fac, fzsc_fac
      USE precon2d, ONLY: ljog_test, lsweep_fast
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(out) :: frzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           frcc, frss, fzcs, fzsc, flcs, flsc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: tempr, tempz
      REAL(rprec) :: t1, t2
      LOGICAL :: lflam
C-----------------------------------------------
      frcc => frzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      fzsc => frzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      flsc => frzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN 
         frss => frzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         fzcs => frzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         flcs => frzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      nsz = ns*nzeta
      lflam = (iequi.ne.2 .or. ljog_test .or. lasym)   

      ALLOCATE (work1(nsz,12), tempr(nsz), tempz(nsz),
     1          stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 tomnsps'

      ioff = LBOUND(frcc,2)
      joff = LBOUND(frcc,3)

      frzl_array = 0

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1

!
!     BEGIN FOURIER TRANSFORM
!
!       FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
!       FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
!       FLmn =      - d(BLmn)/du + d(CLmn)/dv
!
!       NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         IF (iequi.eq.2 .and. lsweep_fast .and. m.eq.1) j2 = 1   !need for fast, symmetric Hessian
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
         l = 0
         DO i = 1, ntheta2
            jll = l+1;  nsl = nsz+l
            l = l + nsz
            tempr(:) = armn(jll:nsl,mparity) 
     1               + xmpq(m,1)*arcon(jll:nsl,mparity)
            tempz(:) = azmn(jll:nsl,mparity) 
     1               + xmpq(m,1)*azcon(jll:nsl,mparity)
            work1(:,1) = work1(:,1) + tempr(:)*cosmui(i,m) 
     1                              + brmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,7) = work1(:,7) + tempz(:)*sinmui(i,m)
     1                              + bzmn(jll:nsl,mparity)*cosmumi(i,m)
            IF (lflam) 
     1      work1(:,11)= work1(:,11)+ blmn(jll:nsl,mparity)*cosmumi(i,m)
 
            IF (.not.lthreed) CYCLE

            work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,3) = work1(:,3) + tempr(:)*sinmui(i,m)
     1                              + brmn(jll:nsl,mparity)*cosmumi(i,m)
            work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
            work1(:,5) = work1(:,5) + tempz(:)*cosmui(i,m)
     1                              + bzmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)

            IF (.not.lflam) CYCLE

            work1(:,9) = work1(:,9) + blmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,1)*cosnv(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,7)*cosnv(k,n)
               IF (lflam)
     1         flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     2                           + work1(jll:nsl,11)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,2)*sinnvn(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,8)*sinnvn(k,n)
               frss(j2:jmax,ni,mj) = frss(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,3)*sinnv(k,n) 
     2                             + work1(j2l:jmaxl,4)*cosnvn(k,n)
               fzcs(j2:jmax,ni,mj) = fzcs(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,5)*sinnv(k,n)
     2                             + work1(j2l:jmaxl,6)*cosnvn(k,n)

               IF (.not.lflam) CYCLE

               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,12)*sinnvn(k,n)
               flcs(jl:ns,ni,mj) = flcs(jl:ns,ni,mj)
     1                           + work1(jll:nsl,9)*sinnv(k,n)
     2                           + work1(jll:nsl,10)*cosnvn(k,n)
            END DO
         END DO
      END DO

      IF (liota) THEN
         ni = 0+ioff;  mj = 0+joff
         DO jl = 2, ns
            t1 = phip(2)*p5*(icurv(jl) + icurv(jl+1))
            t2 = SUM(clmn(jl:nrzt:ns,0)*wint(2:nrzt:ns))
            flsc(jl, ni, mj) = mscale(0)*nscale(0)*
     1                         (-phip(2)*p5*(icurv(jl) + icurv(jl+1))
     2                       + SUM(clmn(jl:nrzt:ns,0)*wint(2:nrzt:ns)))
         END DO
      END IF

!
!     MAKE R,Z(m=1,n=0) SATISFY AVERAGE FORCE BALANCE EXACTLY
!     NOTE: for m=1, FR ~ Z1*(f0 + f2), FZ ~ R1*(f0 - f2), WHERE
!     f0 is the m=0 component of frho, f2 is m=2 component.
!
!     FOR lsweep_fast, NEED ORIGINAL FORCES SO THIS MUST BE SKIPPED
!
      IF (lforbal .and. .not.(lsweep_fast .and. iequi.eq.2) .and.
     1   ((fsqr+fsqz+fsql).lt.1.E-01_dp .and. iter2.gt.1)) THEN
         mj = 1 + joff
         ni = 0 + ioff
         DO jl = 2, ns-1
            work1(jl,1) = frcc_fac(jl)*frcc(jl,ni,mj)
     1                  + fzsc_fac(jl)*fzsc(jl,ni,mj)
            frcc(jl,ni,mj) = rzu_fac(jl)*(equif(jl)/4 + work1(jl,1))
            fzsc(jl,ni,mj) = rru_fac(jl)*(equif(jl)/4 - work1(jl,1))
         END DO

      END IF
            
      DEALLOCATE (work1, tempr, tempz)

      END SUBROUTINE tomnsps
