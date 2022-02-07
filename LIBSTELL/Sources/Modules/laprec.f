      MODULE laprec

      USE stel_kinds, ONLY: rprec
      INTEGER, PARAMETER :: SP=KIND(1.0), DP=KIND(1.0D0)
      LOGICAL, PARAMETER :: ldouble = (rprec.eq.dp)

      INTERFACE bla_axpy

         SUBROUTINE saxpy(N, SA, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(SP), INTENT(IN) :: SA, SX(*)
            REAL(SP), INTENT(INOUT) :: SY(*)
         END SUBROUTINE saxpy

         SUBROUTINE daxpy(N, SA, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(DP), INTENT(IN) :: SA, SX(*)
            REAL(DP), INTENT(INOUT) :: SY(*)
         END SUBROUTINE daxpy

      END INTERFACE

      INTERFACE bla_copy

         SUBROUTINE scopy (N, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(SP), DIMENSION(*), INTENT(IN) :: SX
            REAL(SP), DIMENSION(*), INTENT(OUT) :: SY
         END SUBROUTINE scopy

         SUBROUTINE dcopy (N, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(DP), DIMENSION(*), INTENT(IN) :: SX
            REAL(DP), DIMENSION(*), INTENT(OUT) :: SY
         END SUBROUTINE dcopy

      END INTERFACE

      INTERFACE bla_dot
 
         FUNCTION sdot (n, sx, incx, sy, incy)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: n, incx, incy
            REAL(SP), INTENT(IN) :: sx(*), sy(*)
            REAL(SP) :: sdot
         END FUNCTION sdot

         FUNCTION ddot (n, sx, incx, sy, incy)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: n, incx, incy
            REAL(DP), INTENT(IN) :: sx(*), sy(*)
            REAL(DP) :: ddot
         END FUNCTION ddot

      END INTERFACE

      INTERFACE bla_gemm

         SUBROUTINE sgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     1                     B, LDB, BETA, C, LDC )
            INTEGER, PARAMETER :: SP=KIND(1.0)
            CHARACTER*1         TRANSA, TRANSB
            INTEGER             M, N, K, LDA, LDB, LDC
            REAL(SP) :: ALPHA, BETA
            REAL(SP) :: A( LDA, * ), B( LDB, * ), C( LDC, * )
         END SUBROUTINE sgemm

         SUBROUTINE dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     1                     B, LDB, BETA, C, LDC )
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            CHARACTER*1         TRANSA, TRANSB
            INTEGER             M, N, K, LDA, LDB, LDC
            REAL(DP) :: ALPHA, BETA
            REAL(DP) :: A( LDA, * ), B( LDB, * ), C( LDC, * )
         END SUBROUTINE dgemm

      END INTERFACE

      INTERFACE bla_ger

         SUBROUTINE sger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: M,N,INCX,INCY,LDA
            REAL(SP), INTENT(IN) :: ALPHA, X(*), Y(*)
            REAL(SP), INTENT(INOUT) :: A(LDA,*)
         END SUBROUTINE sger

         SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: M,N,INCX,INCY,LDA
            REAL(DP), INTENT(IN) :: ALPHA, X(*), Y(*)
            REAL(DP), INTENT(INOUT) :: A(LDA,*)
         END SUBROUTINE dger

      END INTERFACE

      INTERFACE bla_iamax

         INTEGER FUNCTION isamax(n, dx, incx)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: incx, n
            REAL(SP), INTENT(IN) :: dx(*)
         END FUNCTION isamax

         INTEGER FUNCTION idamax(n, dx, incx)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: incx, n
            REAL(DP), INTENT(IN) :: dx(*)
         END FUNCTION idamax

      END INTERFACE

      INTERFACE bla_nrm2

         FUNCTION snrm2 (N, SX, INCX)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER N, INCX
            REAL(SP), DIMENSION(N) :: SX
            REAL(SP) :: snrm2
         END FUNCTION snrm2

         FUNCTION dnrm2 (N, SX, INCX)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER N, INCX
            REAL(DP), DIMENSION(N) :: SX
            REAL(DP) :: dnrm2
         END FUNCTION dnrm2

      END INTERFACE

      INTERFACE bla_scal

         SUBROUTINE sscal (N, SA, SX, INCX)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, INCX
            REAL(SP), INTENT(IN) :: SA
            REAL(SP), INTENT(INOUT) :: SX(*)
         END SUBROUTINE sscal

         SUBROUTINE dscal (N, SA, SX, INCX)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, INCX
            REAL(DP), INTENT(IN) :: SA
            REAL(DP), INTENT(INOUT) :: SX(*)
         END SUBROUTINE dscal

      END INTERFACE

      INTERFACE bla_swap

         SUBROUTINE sswap (N, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(SP), INTENT(INOUT) :: SX(*), SY(*)
         END SUBROUTINE sswap

         SUBROUTINE dswap (N, SX, INCX, SY, INCY)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, INCX, INCY
            REAL(DP), INTENT(INOUT) :: SX(*), SY(*)
         END SUBROUTINE dswap

      END INTERFACE

      INTERFACE bla_trsm

         SUBROUTINE strsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1                    A, LDA, B, LDB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: M, N, LDA, LDB
            REAL(SP), INTENT(IN) :: ALPHA
            CHARACTER, INTENT(IN) :: SIDE, UPLO, TRANSA, DIAG
            REAL(SP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(SP), DIMENSION(LDB,*), INTENT(INOUT) :: B
         END SUBROUTINE strsm

         SUBROUTINE dtrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1                    A, LDA, B, LDB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: M, N, LDA, LDB
            REAL(DP), INTENT(IN) :: ALPHA
            CHARACTER, INTENT(IN) :: SIDE, UPLO, TRANSA, DIAG
            REAL(DP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(DP), DIMENSION(LDB,*), INTENT(INOUT) :: B
         END SUBROUTINE dtrsm

      END INTERFACE

      INTERFACE la_getrf

         SUBROUTINE sgetrf (M, N, A, LDA, IPIV, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: M, N, LDA
            INTEGER, INTENT(OUT) :: INFO
            INTEGER, DIMENSION(*), INTENT(OUT) :: IPIV
            REAL(SP), DIMENSION(LDA,*), INTENT(INOUT) :: A
         END SUBROUTINE sgetrf

         SUBROUTINE dgetrf (M, N, A, LDA, IPIV, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: M, N, LDA
            INTEGER, INTENT(OUT) :: INFO
            INTEGER, DIMENSION(*), INTENT(OUT) :: IPIV
            REAL(DP), DIMENSION(LDA,*), INTENT(INOUT) :: A
         END SUBROUTINE dgetrf

      END INTERFACE

      INTERFACE la_getrs

         SUBROUTINE sgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
            INTEGER, INTENT(OUT) :: INFO
            CHARACTER, INTENT(IN) :: TRANS
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(SP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(SP), DIMENSION(LDB,*), INTENT(INOUT) :: B
         END SUBROUTINE sgetrs

         SUBROUTINE dgetrs(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
            INTEGER, INTENT(OUT) :: INFO
            CHARACTER, INTENT(IN) :: TRANS
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(DP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(DP), DIMENSION(LDB,*), INTENT(INOUT) :: B
         END SUBROUTINE dgetrs

         SUBROUTINE sgetrs1(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
            INTEGER, INTENT(OUT) :: INFO
            CHARACTER, INTENT(IN) :: TRANS
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(SP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(SP), DIMENSION(LDB), INTENT(INOUT) :: B
         END SUBROUTINE sgetrs1

         SUBROUTINE dgetrs1(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
            INTEGER, INTENT(OUT) :: INFO
            CHARACTER, INTENT(IN) :: TRANS
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(DP), DIMENSION(LDA,*), INTENT(IN) :: A
            REAL(DP), DIMENSION(LDB), INTENT(INOUT) :: B
         END SUBROUTINE dgetrs1

      END INTERFACE

      INTERFACE la_gesv

         SUBROUTINE sgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER :: INFO, LDA, LDB, N, NRHS
            INTEGER :: IPIV(N)
            REAL(SP) :: A(LDA, N), B(LDB, NRHS)
         END SUBROUTINE sgesv

         SUBROUTINE dgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER :: INFO, LDA, LDB, N, NRHS
            INTEGER :: IPIV(N)
            REAL(DP) :: A(LDA, N), B(LDB, NRHS)
         END SUBROUTINE dgesv

      END INTERFACE

      INTERFACE la_getf2

         SUBROUTINE sgetf2(M, N, A, LDA, IPIV, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER M, N, LDA, INFO
            INTEGER, DIMENSION(*) :: IPIV
            REAL(SP), DIMENSION(LDA,*) :: A
         END SUBROUTINE sgetf2

         SUBROUTINE dgetf2(M, N, A, LDA, IPIV, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER M, N, LDA, INFO
            INTEGER, DIMENSION(*) :: IPIV
            REAL(DP), DIMENSION(LDA,*) :: A
         END SUBROUTINE dgetf2

      END INTERFACE

      INTERFACE laswp

         SUBROUTINE slaswp(N, A, LDA, K1, K2, IPIV, INCX)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: N, LDA, K1, K2, INCX
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(SP), DIMENSION(LDA,*), INTENT(INOUT) :: A
         END SUBROUTINE slaswp

         SUBROUTINE dlaswp(N, A, LDA, K1, K2, IPIV, INCX)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: N, LDA, K1, K2, INCX
            INTEGER, DIMENSION(*), INTENT(IN) :: IPIV
            REAL(DP), DIMENSION(LDA,*), INTENT(INOUT) :: A
         END SUBROUTINE dlaswp

      END INTERFACE

      INTERFACE li_gbfa

         SUBROUTINE sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: ABD
         END SUBROUTINE sgbfa

         SUBROUTINE dgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: ABD
         END SUBROUTINE dgbfa

         SUBROUTINE sgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: ABD
         END SUBROUTINE sgbfa1

         SUBROUTINE dgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: ABD
         END SUBROUTINE dgbfa1

      END INTERFACE

      INTERFACE li_gbsl

         SUBROUTINE sgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: ABD
            REAL(SP), DIMENSION(N) :: B
         END SUBROUTINE sgbsl

         SUBROUTINE sgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: ABD
            REAL(SP), DIMENSION(N) :: B
         END SUBROUTINE sgbsl1

         SUBROUTINE dgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: ABD
            REAL(DP), DIMENSION(N) :: B
         END SUBROUTINE dgbsl

         SUBROUTINE dgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: ABD
            REAL(DP), DIMENSION(N) :: B
         END SUBROUTINE dgbsl1

      END INTERFACE

      INTERFACE li_gefa

         SUBROUTINE sgefa (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: A
         END SUBROUTINE sgefa

         SUBROUTINE dgefa (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: A
         END SUBROUTINE dgefa

         SUBROUTINE sgefa1 (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: A
         END SUBROUTINE sgefa1

         SUBROUTINE dgefa1 (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: A
         END SUBROUTINE dgefa1

      END INTERFACE

      INTERFACE li_gesl

         SUBROUTINE sgesl (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(SP), DIMENSION(LDA,N), INTENT(IN) :: A
            REAL(SP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE sgesl

         SUBROUTINE sgesl1 (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(SP), DIMENSION(*), INTENT(IN) :: A
            REAL(SP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE sgesl1

         SUBROUTINE dgesl (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(DP), DIMENSION(LDA,N), INTENT(IN) :: A
            REAL(DP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE dgesl

         SUBROUTINE dgesl1 (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(DP), DIMENSION(*), INTENT(IN) :: A
            REAL(DP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE dgesl1

      END INTERFACE

      END MODULE laprec
