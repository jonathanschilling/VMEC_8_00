      SUBROUTINE free_mem_funct3d
      USE vmec_main
      USE realspace
      USE vforces
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0
C-----------------------------------------------

      IF (ALLOCATED(armn))
     1   DEALLOCATE (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,
     1   r1, ru, rv, z1, zu, zv, rcon, zcon, ru0, zu0,
     2   rcon0, zcon0, guu, guv, gvv, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error in funct3d'

      IF (ALLOCATED(brv))
     1   DEALLOCATE (brv, bphiv, bzv, bpolvac, bsqvac, stat=istat1)
      IF (istat1 .ne. 0) STOP 'deallocation error in funct3d'

      END SUBROUTINE free_mem_funct3d
