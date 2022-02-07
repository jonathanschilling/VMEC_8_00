      SUBROUTINE totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE vmec_persistent
      USE precon2d, ONLY: lprec2d, rzl_save, loglam_blks
      USE xstuff, ONLY: xc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni, nsz
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
C-----------------------------------------------
!
!     WORK1        Array of inverse transforms in toroidal angle (zeta), for all radial positions
!
!     ONLY WORKING FOR NOW WITH lasym = F; fix later
      IF (iequi.eq.2 .and. .not.loglam_blks .and. .not.lasym) RETURN

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

      nsz = ns*nzeta
      ALLOCATE (work1(nsz,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

      r11 = 0;  ru1 = 0;  rv1 = 0;  rcn1 = 0
      z11 = 0;  zu1 = 0;  zv1 = 0;  zcn1 = 0
      lu1(:,0) = 1;  lu1(:,1) = 0;  lv1 = 0

!
!     EXTRAPOLATION AT JS=1 FOR M=1 MODES
!     NOTE: THE TWO-POINT EXTRAPOLATION SEEMS TO WORK BETTER
!     AND PRODUCE IMPROVED CONVERGENCE AT LARGE NS 
!     HOWEVER, IT CAN NOT BE USED TO COMPUTE THE TRI-DIAG 2D PRECONDITIONER
!
      IF (iequi.eq.2) THEN
         IF (.not.ALLOCATED(rzl_save)) THEN
            ALLOCATE(rzl_save(0:ntor,1:2*ntmax))
            rzl_save(:,1:2*ntmax) = rzl_array(3,:,m1,1:2*ntmax)
     1                            - rzl_array(2,:,m1,1:2*ntmax)
         END IF
         rzl_array(1,:,m1,1:2*ntmax) = rzl_array(2,:,m1,1:2*ntmax) 
     1                               - rzl_save(:,1:2*ntmax)
      ELSE IF (lprec2d) THEN
         rzl_array(1,:,m1,1:2*ntmax) = rzl_array(2,:,m1,1:2*ntmax)
     1                               - rzl_save(:,1:2*ntmax)
      ELSE
         rzl_array(1,:,m1,1:2*ntmax) = 2*rzl_array(2,:,m1,1:2*ntmax)
     1                               -   rzl_array(3,:,m1,1:2*ntmax)
      END IF

      rzl_array(1,:,m1,2*ntmax+1:) = rzl_array(2,:,m1,2*ntmax+1:)

!
!     EXTRAPOLATION OF M=0 MODES FOR LAMBDA (NOTE: Hessian was computed
!     correctly from analytic formulas, so 2-pt extrapolation ok)
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = 2*lmncs(2,:,m0+joff) - lmncs(3,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     1                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,6) = work1(j1l:nsl,6)
     1                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     1                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
               
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)

               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)

               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         l = 0
         DO i = 1, ntheta2
            j1l = l+1;  nsl = nsz+l
            l = l + nsz
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
           
            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,1)*cosmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,1)*sinmum(i,m)
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,1)*cosmux
            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,6)*sinmu(i,m)

            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,6)*cosmum(i,m)
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,6)*sinmux
             
            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,10)*cosmum(i,m)

            IF (.not.lthreed) CYCLE

            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,2)*sinmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,2)*cosmum(i,m)
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,2)*sinmux
            rv1(j1l:nsl,mparity)  = rv1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,3)*cosmu(i,m) 
     1                            + work1(1:nsz,4)*sinmu(i,m)
            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,5)*cosmu(i,m)

            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,5)*sinmum(i,m)
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,5)*cosmux
            zv1(j1l:nsl,mparity)  = zv1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,7)*cosmu(i,m) 
     1                            + work1(1:nsz,8)*sinmu(i,m)

            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,9)*sinmum(i,m)
            lv1(j1l:nsl,mparity)  = lv1(j1l:nsl,mparity)  
     1                            - (work1(1:nsz,11)*cosmu(i,m) 
     1                            +  work1(1:nsz,12)*sinmu(i,m))
         END DO

      END DO

      DEALLOCATE (work1)

      z01(1:ns) = zmnsc(1:ns,n0+ioff,m1+joff)
      r01(1:ns) = rmncc(1:ns,n0+ioff,m1+joff)
      IF (r01(1) .eq. zero) STOP 'r01(0) = 0 in totzsps'
      dkappa = z01(1)/r01(1)

      END SUBROUTINE totzsps

      SUBROUTINE totzsps_hess(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE vmec_persistent
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
      LOGICAL :: logl, logr, logz
C-----------------------------------------------
!
!     SAME AS totzsps, BUT ONLY COMPUTES perturbation (rzl_array is perturbation) TO EXISTING R,Z FOR
!     A PARTICULAR m,n, ntype VALUE
!
!     
!     STORE (FOR loglam_blks) OR RESTORE VARIABLES NOT RECOMPUTED IN TOTZSP DURING HESSIAN JOGS
!
      IF (lasym) RETURN         !NOT IMPLEMENTED YET FOR ASYMMETRIC CASE
      IF (loglam_blks) THEN
         lu_save = lu1;  lv_save = lv1;
         ru_save = ru1;  rv_save = rv1; r1_save = r11; rcon_save = rcn1
         zu_save = zu1;  zv_save = zv1; z1_save = z11; zcon_save = zcn1
         RETURN
      ELSE
         lu1 = lu_save;  lv1 = lv_save
         ru1 = ru_save;  rv1 = rv_save; r11 = r1_save; rcn1 = rcon_save
         zu1 = zu_save;  zv1 = zv_save; z11 = z1_save; zcn1 = zcon_save
      END IF

      logr = ntype_2d .le. ntmax
      logz = .not.logr .and. (ntype_2d .le.2*ntmax)
      logl = ntype_2d .gt. 2*ntmax

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

!
!     EXTRAPOLATION AT JS=1 FOR M=1 MODES
!     NOTE: THE TWO-POINT EXTRAPOLATION SEEMS TO WORK BETTER
!     HOWEVER, IT CAN NOT BE USED WITH THE TRI-DIAG 2D PRECONDITIONER
!
!     FOR THE FAST JOG METHOD (lsweep_fast), MUST MAINTAIN SYMMETRY
!     SO THAT JS=1 JOG IS INDEPENDENT (NOT TIED) TO JS=2 JOG
!     FOR REGULAR JOG METHOD - WHICH DOES NOT ASSUME SYMMETRY - JS=1 AND JS=2 ARE TIED
!
      IF (.not.lsweep_fast) THEN
         rzl_array(1,:,m1,1:2*ntmax) = rzl_array(2,:,m1,1:2*ntmax)   
      END IF

      rzl_array(1,:,m1,2*ntmax+1:) = rzl_array(2,:,m1,2*ntmax+1:)

!
!     ENFORCE CONSTRAINT ON m=1 MODES FOR 3D (CONSISTENT WITH gcz(zcs) = 0 IN RESIDUE)
!
      IF (lthreed) zmncs(:,:,m1+joff) = 0
!
!     EXTRAPOLATION OF M=0 MODES FOR LAMBDA
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
         m = m_2d
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         n = n_2d
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               IF (logr)
     1         work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     2                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logz)
     1         work1(j1l:nsl,6) = work1(j1l:nsl,6)
     2                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logl)
     1         work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     2                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
     
               IF (logr) THEN          
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logz) THEN
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logl) THEN
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

            END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            cosmux
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            sinmux
            END IF
 
            IF (logl)
     1      lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     2            cosmum(i,m)

            IF (.not.lthreed) CYCLE

            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1            cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1            sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1            cosmu(i,m) + work1(:,4)*sinmu(i,m)
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1            sinmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1            cosmux
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1            cosmu(i,m) + work1(:,8)*sinmu(i,m)
            END IF

            IF (logl) THEN
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1            sinmum(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1            cosmu(i,m) + work1(:,12)*sinmu(i,m))
            END IF

         END DO
      
!     PUT THIS IN TOTZSPA_HESS, EVENTUALLY
      IF (m .eq. m1) rzl_array(1,:,m1,:) = 0

      DEALLOCATE (work1)

      END SUBROUTINE totzsps_hess
