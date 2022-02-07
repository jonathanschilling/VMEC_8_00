      MODULE system_mod

      INTERFACE system
         SUBROUTINE vmec_system(cmd, error)
         CHARACTER*(*), INTENT(in) :: cmd
         INTEGER, OPTIONAL :: error
         END SUBROUTINE vmec_system
      END INTERFACE
      
      INTERFACE chdir
         INTEGER FUNCTION vmec_chdir(path)
         CHARACTER*(*), INTENT(in) :: path
         END FUNCTION vmec_chdir
      END INTERFACE
      
      INTERFACE getenv
         SUBROUTINE vmec_getenv(ename, evalue)
         CHARACTER*(*) :: ename, evalue
         END SUBROUTINE vmec_getenv
      END INTERFACE

      INTERFACE putenv
         SUBROUTINE vmec_putenv(ename, evalue, ierror)
         CHARACTER*(*) :: ename, evalue
         INTEGER :: ierror
         END SUBROUTINE vmec_putenv
      END INTERFACE

      INTERFACE PXFFORK
          SUBROUTINE pxffork_g (ipid, ierror)
          INTEGER :: ipid, ierror
          END SUBROUTINE pxffork_g
      END INTERFACE

      INTERFACE PXFGETPID
          SUBROUTINE pxfgetpid_g(ipid, ierror)
          INTEGER :: ipid, ierror
          END SUBROUTINE pxfgetpid_g
      END INTERFACE

      INTERFACE PXFWAIT
          SUBROUTINE pxfwait_g(istat, iretpid, ierror)
          INTEGER :: istat, iretpid, ierror
          END SUBROUTINE pxfwait_g
      END INTERFACE

      END MODULE system_mod
