      SUBROUTINE symforce(ars, brs, crs, azs, bzs, czs, bls, cls, rcs,
     1   zcs, ara, bra, cra, aza, bza, cza, bla, cla, rca, zca)
      USE vmec_main, p5 => cp5
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: ars, brs, crs, azs, bzs, czs,
     2   bls, cls, rcs, zcs
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::
     1   ara, bra, cra, aza, bza, cza, bla, cla, rca, zca
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mpar, ir, i, jk, jka
      REAL(rprec), DIMENSION(ns*nzeta) :: ars_0, brs_0, azs_0, bzs_0,
     1            bls_0, rcs_0, zcs_0, crs_0, czs_0, cls_0
C-----------------------------------------------

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * COS(mu - nv)
!       ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * SIN(mu - nv)
!
!
      DO mpar = 0, 1
         DO i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            IF (i == 1) ir = 1
            DO jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               ara(jk,i,mpar) = p5*(ars(jk,i,mpar)-ars(jka,ir,mpar))
               ars_0(jk)      = p5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
               bra(jk,i,mpar) = p5*(brs(jk,i,mpar)+brs(jka,ir,mpar))
               brs_0(jk)      = p5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
               aza(jk,i,mpar) = p5*(azs(jk,i,mpar)+azs(jka,ir,mpar))
               azs_0(jk)      = p5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
               bza(jk,i,mpar) = p5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar))
               bzs_0(jk)      = p5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
               bla(jk,i,mpar) = p5*(bls(jk,i,mpar)-bls(jka,ir,mpar))
               bls_0(jk)      = p5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
               rca(jk,i,mpar) = p5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar))
               rcs_0(jk)      = p5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
               zca(jk,i,mpar) = p5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar))
               zcs_0(jk)      = p5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
            END DO

            ars(:,i,mpar) = ars_0(:)
            brs(:,i,mpar) = brs_0(:)
            azs(:,i,mpar) = azs_0(:)
            bzs(:,i,mpar) = bzs_0(:)
            bls(:,i,mpar) = bls_0(:)
            rcs(:,i,mpar) = rcs_0(:)
            zcs(:,i,mpar) = zcs_0(:)

            IF (lthreed) THEN
               DO jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  cra(jk,i,mpar) = p5*(crs(jk,i,mpar)+crs(jka,ir,mpar))
                  crs_0(jk)      = p5*(crs(jk,i,mpar)-crs(jka,ir,mpar))
                  cza(jk,i,mpar) = p5*(czs(jk,i,mpar)-czs(jka,ir,mpar))
                  czs_0(jk)      = p5*(czs(jk,i,mpar)+czs(jka,ir,mpar))
                  cla(jk,i,mpar) = p5*(cls(jk,i,mpar)-cls(jka,ir,mpar))
                  cls_0(jk)      = p5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
               END DO

               crs(:,i,mpar) = crs_0(:)
               czs(:,i,mpar) = czs_0(:)
               cls(:,i,mpar) = cls_0(:)
            ENDIF

         END DO
      END DO

      END SUBROUTINE symforce
