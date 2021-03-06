      SUBROUTINE allocate_funct3d
      USE vmec_main
      USE realspace
      USE vforces
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1, ndim, ndim2
C-----------------------------------------------

      ndim  = 1+nrzt
      ndim2 = 2*ndim

      CALL free_mem_funct3d

      ALLOCATE( armn(ndim2), azmn(ndim2), brmn(ndim2), bzmn(ndim2),
     1   crmn(ndim2), czmn(ndim2), blmn(ndim2), clmn(ndim2),
     2   r1(nrzt,0:1), ru(nrzt,0:1), rv(nrzt,0:1),
     3   z1(nrzt,0:1), zu(nrzt,0:1), zv(nrzt,0:1),
     4   rcon(nrzt,0:1), zcon(nrzt,0:1), ru0(ndim), zu0(ndim),
     6   rcon0(ndim), zcon0(ndim), guu(ndim), guv(ndim), gvv(ndim),
     7   stat=istat1 )
      IF (istat1.ne.0) STOP 'allocation error #1 in funct3d'

      IF (lfreeb) THEN
      ALLOCATE (brv(nznt), bphiv(nznt), bzv(nznt),
     1   bpolvac(nznt), bsqvac(nznt), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #2 in funct3d'
      END IF
!
!     Pointer alias assignments
!     NOTE: In FORCES, X_e(nrzt+1) overlaps X_o(1), which should never be used...
!
      armn_e => armn(:ndim)
      armn_o => armn(ndim:)
      armn(:ndim2) = zero
      brmn_e => brmn(:ndim)
      brmn_o => brmn(ndim:)
      brmn(:ndim2) = zero
      azmn_e => azmn(:ndim)
      azmn_o => azmn(ndim:)
      azmn(:ndim2) = zero
      bzmn_e => bzmn(:ndim)
      bzmn_o => bzmn(ndim:)
      bzmn(:ndim2) = zero
      crmn_e => crmn(:ndim)
      crmn_o => crmn(ndim:)
      crmn(:ndim2) = zero
      czmn_e => czmn(:ndim)
      czmn_o => czmn(ndim:)
      czmn(:ndim2) = zero
      blmn_e => blmn(:ndim)
      blmn_o => blmn(ndim:)
      blmn(:ndim2) = zero
      clmn_e => clmn(:ndim)
      clmn_o => clmn(ndim:)
      clmn(:ndim2) = zero
      rcon0(:ndim) = zero
      zcon0(:ndim) = zero

      END SUBROUTINE allocate_funct3d
