      SUBROUTINE jxbforce(bsupu, bsupv, bsubu, bsubv, bsubs, bsubsu,
     1   bsubsv, itheta, brho, gsqrt, izeta, bsq)
      USE safe_open_mod
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, signgs, mnyq, nnyq
      USE realspace
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1  bsupu, bsupv, bsq, gsqrt
      REAL(rprec), DIMENSION(ns,nznt,0:1), INTENT(inout) ::
     1  bsubu, bsubv
      REAL(rprec), DIMENSION(ns,nznt), INTENT(inout) :: bsubs
      REAL(rprec), DIMENSION(ns,nznt), INTENT(out) ::
     1  itheta, brho, izeta
      REAL(rprec), DIMENSION(ns,nznt,0:1) :: bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      LOGICAL, PARAMETER :: lbsubs = .true.       !!True  to USE NEW bsubs calculation (from mag. diff. eq.)
                                                  !!False to USE OLD bsubs calculation (from metrics)
      LOGICAL, PARAMETER :: lprint = .false.      !!Prints out bsubs spectrum to fort.33
      REAL(rprec), PARAMETER :: two=2, p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER lk, lz, lt, k, m, js, j, n, injxbout, mparity
      INTEGER :: njxbout = jxbout0, nmin
      INTEGER, PARAMETER :: ns_skip = 1, nu_skip = 1, nv_skip = 1
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bdotj, bsubuv, bsubvu, brhomn
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmn
      REAL(rprec), DIMENSION(:), ALLOCATABLE     :: jperpu, jperpv, 
     2    sqgb2, sqrtg, jp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: bsubua, bsubva
      REAL(rprec) ::
     1    bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,
     1    bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,
     2    dnorm1, tsini1, tsini2, tcosi1, tcosi2, tcosm1, tcosm2,
     3    tcosn1, tcosn2, tsinm1, tsinm2, tcos1, tcos2, tsin1, tsin2,
     4    tsinn1, tsinn2, avforce, pprime, amaxfor, aminfor, tjnorm,
     5    ovp, pnorm, brho00(ns)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a
      CHARACTER :: jxbout_file*100
C-----------------------------------------------

      jxbout_file = 'jxbout.'//input_extension
      CALL safe_open(njxbout, injxbout, jxbout_file, 'replace',
     1    'formatted')
      IF (injxbout .ne. 0) THEN
         PRINT *,' Error opening JXBOUT file in jxbforce'
         RETURN
      END IF

      WRITE (njxbout,6) (ns1-1)/ns_skip, ntheta2/nu_skip, nzeta/nv_skip,
     1    mpol, ntor
 6    FORMAT(/,' Radial surfaces = ',i3, ' Poloidal grid points = ',i3,
     1         ' Toroidal grid points = ',i3,/,
     2         ' Poloidal modes = ',i3,' Toroidal modes = ', i3)

!
!     PROGRAM FOR COMPUTING LOCAL JXB = grad-p FORCE BALANCE
!
!     Compute u (=theta), v (=zeta) derivatives of B sub s
!
      WRITE (njxbout, 5)
 5    FORMAT(' LEGEND:',/,
     1  " U = VMEC poloidal angle, V = VMEC (geometric) toroidal angle"/
     2  " SQRT(g') = SQRT(g-VMEC) / V': Jacobian based on VOLUME",/,
     3  " V' = dV/ds: differential volume element",/,
     4  " Es = SQRT(g') [grad(U) X grad(V)] : covariant radial",
     4  " unit vector based on volume",/,
     5  " BSUP{U,V} = B DOT GRAD{U,V}:  contravariant components of B"/,
     6  " JSUP{U,V} = SQRT(g') J DOT GRAD{U,V}",/,
     7  " J X B = Es DOT [J X B]: covariant component of J X B force",/,
     8  " J * B = J DOT B * SQRT(g')",/,
     9  " p' = dp/dV: pressure gradient based on volume derivative",//)

      lz = nzeta*ntheta2
      ALLOCATE (bdotj(ns,nznt), bsubuv(ns,nznt),
     1          bsubvu(ns,nznt), jperpu(nznt), jperpv(nznt), 
     2          sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq),jp2(nznt),
     3          bsubua(nznt,0:1), bsubva(nznt,0:1), jxb(nznt), 
     4          jxb2(nznt), bsupu1(nznt), bsupv1(nznt), bsubu1(nznt), 
     5          bsubv1(nznt), bsubsmn(ns,0:mnyq,-nnyq:nnyq),
     6          bsubs_s(lz,0:1), bsubs_a(lz,0:1), sqrtg(nznt),
     7          bsubu_s(lz,0:1), bsubu_a(lz,0:1),
     8          bsubv_s(lz,0:1), bsubv_a(lz,0:1), stat = j)
      IF (j .ne. 0) STOP 'allocation error in jxbforce'


!
!     NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
!
      bsubsu = 0; bsubsv = 0; bsubuv = 0; bsubvu = 0; bdotj  = 0

      radial: DO js = 1, ns
!
!     Put bsubs on full mesh
!
         IF (js.gt.1 .and. js.lt.ns) THEN     
            bsubs(js,:) = p5*(bsubs(js,:) + bsubs(js+1,:))
         END IF

         bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
         bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
         bsubua = 0;   bsubva = 0

         IF (lasym)  THEN
            CALL fsym_fft (bsubs(js,:), bsubu(js,:,:), bsubv(js,:,:),
     1             bsubs_s, bsubu_s, bsubv_s, bsubs_a, bsubu_a, bsubv_a)
         ELSE
            bsubs_s(:,0) = bsubs(js,:); bsubu_s = bsubu(js,:,:)
            bsubv_s = bsubv(js,:,:)
         END IF

         DO m = 0, mnyq
            mparity = MOD(m, 2)
            DO n = 0, nnyq
               dnorm1 = one/r0scale**2
               IF (m.eq.mnyq) dnorm1 = p5*dnorm1
               IF (n.eq.nnyq .and. n.ne.0) dnorm1 = p5*dnorm1
               bsubsmn1 = 0;  bsubsmn2 = 0
               IF (lasym) bsubsmn3 = 0;  bsubsmn4 = 0
               IF (m.gt.mpol1 .or. n.gt.ntor) GOTO 222
               bsubumn1 = 0;  bsubumn2 = 0;  bsubvmn1 = 0;  bsubvmn2 = 0
               IF (lasym)
     1         bsubumn3 = 0;  bsubumn4 = 0;  bsubvmn3 = 0;  bsubvmn4 = 0

               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1
                     tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1
                     bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk,0)
                     bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk,0)
                     bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity)
                     bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity)
                     bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity)
                     bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity)

                     IF (lasym) THEN
                     bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk,0)
                     bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk,0)
                     bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity)
                     bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity)
                     bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity)
                     bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity)
                     END IF

                     lk = lk + nzeta
                  END DO
               END DO

!
!              Compute on half u grid (must add symmetric, antisymmetric parts for lasym)
! 
               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tcos1 = cosmu(j,m)*cosnv(k,n)
                     tcos2 = sinmu(j,m)*sinnv(k,n)
                     bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 +
     1                  tcos2*bsubumn2
                     bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 +
     1                  tcos2*bsubvmn2

                     tcosm1 = cosmum(j,m)*cosnv(k,n)
                     tcosm2 = sinmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     bsubvu(js,lk) = bsubvu(js,lk) + (
     1                               sinmum(j,m)*cosnv(k,n)*bsubvmn1 +
     2                               cosmum(j,m)*sinnv(k,n)*bsubvmn2)
                     bsubuv(js,lk) = bsubuv(js,lk) + (
     1                               cosmu(j,m)*sinnvn(k,n)*bsubumn1 +
     2                               sinmu(j,m)*cosnvn(k,n)*bsubumn2)

                     IF (lasym) THEN
                     tsin1 = sinmui(j,m)*cosnv(k,n)
                     tsin2 = cosmui(j,m)*sinnv(k,n)
                     bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 +
     1                  tsin2*bsubumn4
                     bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 +
     1                  tsin2*bsubvmn4

                     tsinm1 = sinmum(j,m)*cosnv(k,n)
                     tsinm2 = cosmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                   tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                   tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                     END IF
                     lk = lk + nzeta
                  END DO
               END DO

 222           CONTINUE

!
!        FIX MULTIPLIER FOR M != 0, N != 0 (p5)
!
               IF (m .eq. 0) THEN
                  bsubsmn(js,m,n) =-bsubsmn2
                  IF (n.ne.0) bsubsmn(js,m,-n)= 0
               ELSE IF (n .eq. 0) THEN
                  bsubsmn(js,m,n) = bsubsmn1
               ELSE
                  bsubsmn(js,m,n)  = p5*(bsubsmn1 - bsubsmn2)
                  bsubsmn(js,m,-n) = p5*(bsubsmn1 + bsubsmn2)
               END IF

            END DO
         END DO

         IF (lasym) CALL fsym_invfft (bsubua, bsubva, 1)
         bsubu(js,:,0) = bsubua(:,0)
         bsubv(js,:,0) = bsubva(:,0)

      END DO radial

      IF (lasym) CALL fsym_invfft (bsubsu, bsubsv, ns)

      IF (.not.lbsubs) GOTO 1500          !!SKIPS Bsubs Correction - uses Bsubs from metric elements

!
!     Compute corrected B-sub-s (impacts currents)
!
      correct_bsubs: DO js = 2, ns-1
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
         brho(js,:) = ohs*
     1   ( bsupu1(:)*(bsubu(js+1,:,0) - bsubu(js,:,0))
     2   + bsupv1(:)*(bsubv(js+1,:,0) - bsubv(js,:,0)))
     3   + (pres(js+1) - pres(js))*ohs*jxb(:)
         brho00(js) = SUM(brho(js,:)*wint(js:nrzt:ns))
         brho(js,:) = brho(js,:) - signgs*jxb(:)*brho00(js)/
     1      (p5*(vp(js) + vp(js+1)))

         jxb(:) = brho(js,:)
         CALL getbrho (brhomn, jxb, bsupu1, bsupv1, mnyq, nnyq)
         IF (lprint) THEN
            WRITE (33, *) ' JS = ', js
            WRITE (33, *) '  M    N       BSUBS(old)        BSUBS(new)'
            DO m = 0, mnyq
               nmin = -nnyq
               IF (m .eq. 0) nmin = 0
               DO n = nmin, nnyq
                  WRITE(33,1223) m, n, bsubsmn(js,m,n), brhomn(m,n)
               END DO
            END DO
         END IF
 1223    FORMAT (i4,1x,i4,2(6x,1p,e12.3))

!
!        Recompute bsubsu,v now using corrected bsubs
!        Store old values (itheta,izeta) for checking force balance later
!
         itheta(js,:) = bsubsu(js,:,0);  izeta (js,:) = bsubsv(js,:,0)
         bsubs(js,:) = 0;   bsubsu(js,:,:) = 0;   bsubsv(js,:,:) = 0

         DO m = 0, mnyq
            DO n = 0, nnyq
               IF (n .eq. 0) THEN
                  bsubsmn1 = brhomn(m,n)
                  bsubsmn2 = 0
               ELSE
                  bsubsmn1 = brhomn(m,n) + brhomn(m,-n)
                  bsubsmn2 =-brhomn(m,n) + brhomn(m,-n)
               END IF
               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tsin1 = sinmu(j,m)*cosnv(k,n)
                     tsin2 = cosmu(j,m)*sinnv(k,n)
                     bsubs(js,lk) = bsubs(js,lk) + tsin1*bsubsmn1
     1                                           + tsin2*bsubsmn2
                     tcosm1 = cosmum(j,m)*cosnv(k,n)
                     tcosm2 = sinmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     lk = lk + nzeta
                  END DO
               END DO
            END DO
         END DO

      END DO correct_bsubs

      IF (lasym) CALL fsym_invfft (bsubsu, bsubsv, ns)

!
!     CHECK FORCE BALANCE: SQRT(g)*(bsupu*bsubsu + bsupv*bsubsv) = brho
!
!     check_fb: DO js = 2, ns-1
!        bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
!    1    + bsupu(js+1,:)*gsqrt(js+1,:))
!        bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
!    1    + bsupv(js+1,:)*gsqrt(js+1,:))
!        jp2(:) = bsupu1(:)*bsubsu(js,:,0) + bsupv1(:)*bsubsv(js,:,0)
!        jxb(:) = bsupu1(:)*itheta(js,:) + bsupv1(:)*izeta(js,:)

!        PRINT *,'JS = ',js
!        DO lk = 1, nznt
!           WRITE(*,1224) lk, brho(js,lk), jxb(lk), jp2(lk)
!        END DO

!        PAUSE

!     END DO check_fb

 1224 FORMAT ('lk = ',i4,' brho(rhs) = ', 1p,e12.4,
     1   ' B dot grad Bs(old) = ', 1p,e12.4, ' B dot grad Bs(new) = ',
     2   1p,e12.4)

 1500 CONTINUE

      DEALLOCATE (bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, 
     1            bsubv_a, stat=lk)

!
!     Compute end point values for bsubs
!
      bsubs(1,:)  = 2*bsubs(2,:)  - bsubs(3,:)
      bsubs(ns,:) = 2*bsubs(ns,:) - bsubs(ns-1,:)
!
!     Now compute currents on the FULL radial mesh
!     Here:
!
!     Itheta = sqrt(g) * Jsupu
!     Izeta  = sqrt(g) * Jsupv
!     Jsupx  = J dot grad(x)                          x=(u,v)
!     jxb    = (J X B) dot (grad-u X grad-v) sqrt(g)  
!     bdotj  = sqrt(g)*J dot B
!     jperpx = (B X gradp) dot grad(x) / |B|**2       x=(u,v)
!     sqgb2  = sqrt(g)*|B|**2
!     sqrtg  = sqrt(g)
!     pprime = dp/dV
!
!     jp2   == |j-perp|**2 = jperpu**2 * guu + 2*jperpu*jperpv*guv + jperv**2 * gvv
!     This was compared to the alternative expression (agreed very well):
!     |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
!
!     Note: Multiply currents, pressure by 1/mu0 to get in mks units!
!
      DO js = 2, ns1
         ovp = two/(vp(js+1) + vp(js))
         tjnorm = ovp*signgs
         sqgb2(:nznt) = (gsqrt(js+1,:nznt)*(bsq(js+1,:nznt)- pres(js+1))
     1                +  gsqrt(js,:nznt)  *(bsq(js,:nznt) - pres(js)))
         IF (ANY(sqgb2(:nznt)*signgs .le. zero)) 
     1       STOP ' SQGB2 <= 0 in JXBFORCE'
         pprime = ohs*(pres(js+1)-pres(js))/mu0              !dp/ds here
         jperpu(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))*
     1                       pprime/sqgb2
         jperpv(:nznt) =-p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))*
     1                       pprime/sqgb2
         jp2(:nznt)=p5*(jperpu**2*(guu(js+1:nrzt:ns) + guu(js:nrzt:ns))
     1          + 2*jperpu*jperpv*(guv(js+1:nrzt:ns) + guv(js:nrzt:ns))
     2          +       jperpv**2*(gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))
         itheta(js,:nznt) =  bsubsv(js,:nznt,0) - ohs*
     1                      (bsubv(js+1,:nznt,0) - bsubv(js,:nznt,0))
         izeta(js,:nznt)  = -bsubsu(js,:nznt,0) + ohs*
     1                      (bsubu(js+1,:nznt,0) - bsubu(js,:nznt,0))
         itheta(js,:nznt) = itheta(js,:nznt)/mu0
         izeta(js,:nznt)  = izeta(js,:nznt)/mu0
         sqrtg(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:nznt) = p5*(bsupu(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupu(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)
         bsupv1(:nznt) = p5*(bsupv(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupv(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)
         bsubu1(:nznt) = p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))
         bsubv1(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))
         jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt)
     1              -      izeta (js,:nznt) * bsupu1(:nznt))
         bdotj(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) +
     1                     izeta (js,:nznt) * bsubv1(:nznt)
         pprime = ovp*pprime
         pnorm = one/(ABS(pprime) + EPSILON(pprime))
         amaxfor = MAX(MAXVAL(jxb(1:nznt)-pprime)*pnorm, zero)
         aminfor = MIN(MINVAL(jxb(1:nznt)-pprime)*pnorm, zero)
         avforce = SUM(wint(2:nrzt:ns)*(jxb(:nznt) - pprime))


!        Compute <j dot B>, <B sup v> = signgs*phip
!        jpar2 = <j||**2>, jperp2 = <j-perp**2>,  with <...> = flux surface average

         jdotb(js) = tjnorm*SUM(bdotj(js,:nznt)*wint(2:nrzt:ns))
         bdotb(js) = tjnorm*SUM(sqgb2(:nznt)*wint(2:nrzt:ns))
         bdotgradv(js) = p5*(phip(js) + phip(js+1))*tjnorm
         jpar2(js) = tjnorm*SUM(bdotj(js,:nznt)**2*wint(2:nrzt:ns)
     1                         /sqgb2(:nznt))
         jperp2(js)= tjnorm*SUM(jp2(:nznt)*wint(2:nrzt:ns)*sqrtg(:nznt))

         IF (MOD(js,ns_skip) .eq. 0) THEN
            amaxfor = MIN(amaxfor,9.999_dp)
            aminfor = MAX(aminfor,-9.999_dp)
            WRITE (njxbout, 200) phi(js)/phi(ns), avforce, jdotb(js),
     1         bdotgradv(js), 100.0_dp*amaxfor, 100.0_dp*aminfor
            WRITE (njxbout, 90)
            DO lz = 1, nzeta, nv_skip
               WRITE (njxbout, 100) 360._dp*(lz-1)/nzeta, lz
               DO lt = 1, ntheta2, nu_skip
                  lk = lz + nzeta*(lt - 1)
                  WRITE (njxbout, 110) lt, ovp*itheta(js,lk),
     1              ovp*izeta(js,lk), ovp*(bsubuv(js,lk) -
     2              bsubvu(js,lk))/mu0, bsupu1(lk), bsupv1(lk),
     3              jxb(lk), pprime,
     4        (jxb(lk) - pprime), ovp*bdotj(js,lk), bsubu(js,lk,0),
     5              bsubv(js,lk,0), bsubs(js,lk)
               END DO
            END DO
         ENDIF
      END DO

      CLOSE (njxbout)

      izeta(1,:nznt) = two*izeta(2,:nznt) - izeta(3,:nznt)           !!For print out in wrout
      izeta(ns,:nznt)= two*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)     !!For print out in wrout
      jdotb(1) = two*jdotb(2) - jdotb(3)
      jdotb(ns) = two*jdotb(ns-1) - jdotb(ns-2)
      bdotb(1) = two*bdotb(3) - bdotb(2)
      bdotb(ns) = two*bdotb(ns-1) - bdotb(ns-2)
      bdotgradv(1) = two*bdotgradv(2) - bdotgradv(3)
      bdotgradv(ns) = two*bdotgradv(ns-1) - bdotgradv(ns-2)
      jpar2(1)   = 0; jpar2(ns)  = 0; jperp2(1)  = 0; jperp2(ns) = 0

      DEALLOCATE (jperpu, jperpv, sqgb2, sqrtg, jp2, brhomn, bsubsmn, 
     1    bsubua, bsubva, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1)

!
!     COMPUTE MERCIER CRITERION
!
      bdotj = mu0*bdotj
      CALL Mercier(gsqrt,bsq,bdotj,iotas,wint,r1,ru,rv,zu,zv,bsubu,
     1             vp,phips,pres,ns,nznt)


   90 FORMAT(/"   LU      JSUPU      JSUPV      JSUPS      BSUPU",
     1   "      BSUPV      J X B       p'    J X B - p'     J * B",
     2   "      BSUBU      BSUBV      BSUBS   "/)
  100 FORMAT( " TOROIDAL ANGLE (PER PERIOD) = ", f8.3," DEGREES",
     1        " (PLANE #", i3,")")
  110 FORMAT(i5,1p,12e11.3)
  200 FORMAT(/" TOROIDAL FLUX = ",1p,e12.3,3x,"<J X B - p'> = ",
     1   e12.3,3x,"<J DOT B> = ",e12.3,3x,
     2   "<B DOT GRAD(V)> = ",e12.3,/,
     3   " MAXIMUM FORCE DEVIATIONS (RELATIVE TO p'): ",sp,0p,f7.2,"%",
     4     3x,f7.2,"%")

      DEALLOCATE (bdotj, bsubuv, bsubvu, stat = j)

      END SUBROUTINE jxbforce
