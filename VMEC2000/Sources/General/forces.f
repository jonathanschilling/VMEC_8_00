      SUBROUTINE forces
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces
      USE vsvd, ONLY: torflux_edge => torflux
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, lk, ndim
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    bsqr, gvvs, guvs, guus
      REAL(rprec), DIMENSION(:), POINTER :: gcon, lv_e, lu_e, lu_o
      REAL(rprec) :: rcon1, zcon1, dphids
C-----------------------------------------------
      ndim = 1+nrzt

!
!     POINTER ALIASES
!
      gcon => z1(:,0)
      lv_e => crmn_e; lu_e => czmn_e; lu_o => czmn_o

      ALLOCATE (bsqr(ndim), gvvs(ndim), guvs(ndim), guus(ndim),
     1           stat=l)
      IF (l .ne. 0) STOP 'Allocation error in VMEC FORCES'

!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT .5*(RU-odd*Zodd - ZU-odd*Rodd) TERM
!           IN JACOBIAN
!
!     LV :  VARIATION OF R-TERM IN JACOBIAN
!
!     GVV:  VARIATION OF R**2-TERM AND Rv**2,Zv**2 IN gvv
!
!     GUU, GUV: VARIATION OF Ru, Rv, Zu, Zv IN guu, guv
!
      dphids = p25/torflux_edge

      lu_e(1:ndim:ns) = 0; lv_e(1:ndim:ns) = 0
      guu(1:ndim:ns)  = 0; guv(1:ndim:ns)  = 0; gvv(1:ndim:ns) = 0
      guus = guu*shalf;    guvs = guv*shalf;    gvvs = gvv*shalf

CDIR$ IVDEP
      DO l = 1, ndim
         armn_e(l)  = ohs*armn_e(l) * lu_e(l)
         azmn_e(l)  =-ohs*azmn_e(l) * lu_e(l)
         brmn_e(l)  = brmn_e(l) * lu_e(l)
         bzmn_e(l)  =-bzmn_e(l) * lu_e(l)
         bsqr(l)    = phip(l)*lu_e(l)/shalf(l)
      END DO

CDIR$ IVDEP
      DO l = 1, ndim
         armn_o(l)  = armn_e(l) *shalf(l)
         azmn_o(l)  = azmn_e(l) *shalf(l)
         brmn_o(l)  = brmn_e(l) *shalf(l)
         bzmn_o(l)  = bzmn_e(l) *shalf(l)
      END DO
!
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
CDIR$ IVDEP
      DO l = 1, nrzt
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = dphids*(bsqr(l) + bsqr(l+1))
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
      END DO

CDIR$ IVDEP
      DO l = 1, nrzt
         armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
     1             - gvv(l)*r1(l,0)
         azmn_e(l) = azmn_e(l+1) - azmn_e(l)
         brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
         bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
      END DO

CDIR$ IVDEP
      DO l = 1, nrzt
        armn_e(l) = armn_e(l) - gvvs(l)*r1(l,1)
        brmn_e(l) = brmn_e(l) + bsqr(l)*z1(l,1)
     1            - guus(l)*ru(l,1) - guu(l)*ru(l,0)
        bzmn_e(l) = bzmn_e(l) - bsqr(l)*r1(l,1)
     1            - guus(l)*zu(l,1) - guu(l)*zu(l,0)
      END DO

      lv_e(1:nrzt+1) = lv_e(1:nrzt+1)*shalf(1:nrzt+1)
CDIR$ IVDEP
      DO l = 1, nrzt
         armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l)
     1             + p5*(lv_e(l) + lv_e(l+1))
         azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l)
         brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
         bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))
      END DO

      guu(1:nrzt)  = guu(1:nrzt) * sqrts(1:nrzt)**2
      bsqr(1:nrzt) = gvv(1:nrzt) * sqrts(1:nrzt)**2
CDIR$ IVDEP
      DO l = 1, nrzt
         lu_o(l)   = dphids*(lu_e(l)*phip(l) + lu_e(l+1)*phip(l+1))
         armn_o(l) = armn_o(l) - (zu(l,1)*lu_o(l)
     1             + bsqr(l)*r1(l,1) + gvvs(l)*r1(l,0))
         azmn_o(l) = azmn_o(l) + ru(l,1)*lu_o(l)
         brmn_o(l) = brmn_o(l) + z1(l,1)*lu_o(l)
     1             -(guu(l)*ru(l,1) + guus(l)*ru(l,0))
         bzmn_o(l) = bzmn_o(l) - (r1(l,1)*lu_o(l)
     1             + guu(l)*zu(l,1) + guus(l)*zu(l,0))
      END DO

      IF (lthreed) THEN
CDIR$ IVDEP
         DO l = 1, nrzt
            guv(l)  = p5*(guv(l) + guv(l+1))
            guvs(l) = p5*(guvs(l) + guvs(l+1))
            brmn_e(l) = brmn_e(l) - (guv(l)*rv(l,0) + guvs(l)*rv(l,1))
            bzmn_e(l) = bzmn_e(l) - (guv(l)*zv(l,0) + guvs(l)*zv(l,1))
            crmn_e(l) = guv(l) *ru(l,0) + gvv(l) *rv(l,0)
     1                + gvvs(l)*rv(l,1) + guvs(l)*ru(l,1)
            czmn_e(l) = guv(l) *zu(l,0) + gvv(l) *zv(l,0)
     1                + gvvs(l)*zv(l,1) + guvs(l)*zu(l,1)
         END DO

CDIR$ IVDEP
         DO l = 1, nrzt
            guv(l) = guv(l) *sqrts(l)*sqrts(l)
            brmn_o(l) = brmn_o(l) - (guvs(l)*rv(l,0) + guv(l)*rv(l,1))
            bzmn_o(l) = bzmn_o(l) - (guvs(l)*zv(l,0) + guv(l)*zv(l,1))
            crmn_o(l) = guvs(l)*ru(l,0) + gvvs(l)*rv(l,0)
     1                + bsqr(l)*rv(l,1) + guv(l) *ru(l,1)
            czmn_o(l) = guvs(l)*zu(l,0) + gvvs(l)*zv(l,0)
     1                + bsqr(l)*zv(l,1) + guv(l) *zu(l,1)
         END DO
      ENDIF
!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      IF (ivac .ge. 1) THEN
         lk = 0
CDIR$ IVDEP
         DO l = ns,nrzt,ns
            lk = lk+1
            armn_e(l) = armn_e(l) + zu0(l)*rbsq(lk)
            armn_o(l) = armn_o(l) + zu0(l)*rbsq(lk)
            azmn_e(l) = azmn_e(l) - ru0(l)*rbsq(lk)
            azmn_o(l) = azmn_o(l) - ru0(l)*rbsq(lk)
         END DO
         fz00_edge = SUM(wint(ns:nrzt:ns)*ru0(ns:nrzt:ns)*rbsq(1:nznt))
      ENDIF

 100  CONTINUE

      DEALLOCATE (bsqr, gvvs, guvs, guus, stat=l)
!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
CDIR$ IVDEP
      DO l = 1,nrzt
         rcon1   = (rcon(l,0) - rcon0(l)) * gcon(l)
         zcon1   = (zcon(l,0) - zcon0(l)) * gcon(l)
         brmn_e(l) = brmn_e(l) + rcon1
         bzmn_e(l) = bzmn_e(l) + zcon1
         brmn_o(l) = brmn_o(l)+ rcon1*sqrts(l)
         bzmn_o(l) = bzmn_o(l)+ zcon1*sqrts(l)
         rcon(l,0) =  ru0(l) * gcon(l)
         zcon(l,0) =  zu0(l) * gcon(l)
         rcon(l,1) = rcon(l,0) * sqrts(l)
         zcon(l,1) = zcon(l,0) * sqrts(l)
      END DO

      END SUBROUTINE forces
