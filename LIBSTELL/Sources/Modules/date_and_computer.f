      MODULE date_and_computer
      CHARACTER*(3), DIMENSION(12), PARAMETER :: months =
     1  ( / 'Jan','Feb','Mar','Apr','May','Jun',
     2      'Jul','Aug','Sep','Oct','Nov','Dec' / )
c !DEC$ IF DEFINED (RISC)
c       CHARACTER*(*), PARAMETER :: computer =
c      1' IBM RISC-6000 WORKSTATION. '
c !DEC$ ELSEIF DEFINED (IRIX64)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' ORIGIN 2000. '
c !DEC$ ELSEIF DEFINED (IRIX)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' PPPL CCF. '
c !DEC$ ELSEIF DEFINED (HPUX)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' HP-UX WORKSTATION. '
c !DEC$ ELSEIF DEFINED (OSF1)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' DECSTATION AXP. '
c !DEC$ ELSEIF DEFINED (SUN)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' SUN WORKSTATION. '
c !DEC$ ELSEIF DEFINED (VAX)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' DEC VAX. '
c !DEC$ ELSEIF DEFINED (WIN32)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' WINDOWS NT SYSTEM. '
c !DEC$ ELSEIF DEFINED (LINUX)
      CHARACTER*(*), PARAMETER :: computer =
     1  ' LINUX SYSTEM. '
c !DEC$ ELSEIF DEFINED (AXPVMS)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' AXP/VMS. '
c !DEC$ ELSEIF DEFINED (CRAY)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' CRAY SUPER-COMPUTER. '
c !DEC$ ELSEIF DEFINED (SX5)
c       CHARACTER*(*), PARAMETER :: computer =
c      1  ' SX-5 SUPER-COMPUTER. '
c !DEC$ ELSE
c       CHARACTER*(*), PARAMETER :: computer =
c      1  '  '
c !DEC$ ENDIF
      END MODULE date_and_computer
