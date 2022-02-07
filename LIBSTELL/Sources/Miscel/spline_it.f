       SUBROUTINE spline_it(ndata, xdata, ydata, npts, x, y, i_full)
       USE stel_kinds
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ndata, npts, i_full
       REAL(rprec), DIMENSION(ndata), INTENT(IN):: xdata, ydata
       REAL(rprec), DIMENSION(npts), INTENT(OUT):: x, y
       REAL(rprec), DIMENSION(npts) :: dy
       REAL(rprec), DIMENSION(:), ALLOCATABLE:: xfull, yfull,
     1   dyfull, wk
       INTEGER:: iwk, ierr
       LOGICAL:: spline

       ALLOCATE(xfull(ndata+i_full), yfull(ndata+i_full),
     1     dyfull(ndata+i_full), wk(2*(ndata+i_full)) )

       IF (i_full .eq. 1) THEN                                            !! I_FULL= 1, data on HALF_mesh
         xfull(1) = -1
         xfull(ndata+1) = 1
         xfull(2:ndata) = 0.5_dp*(xdata(1:ndata-1)+xdata(2:ndata))
         yfull(1) = ydata(1) + (xdata(1)+1.)*
     1    (ydata(2)-ydata(1))/(xdata(2)-xdata(1))
         yfull(2:ndata) = 0.5_dp*(ydata(1:ndata-1)+ydata(2:ndata))
         yfull(ndata+1) = ydata(ndata) + (1 - xdata(ndata))*
     1    (ydata(ndata-1)-ydata(ndata))/(xdata(ndata-1)-xdata(ndata))
       ELSE                                                              !!I_FULL =0, data on FULL mesh;
         xfull(1:ndata) = xdata(1:ndata)
         yfull(1:ndata) = ydata(1:ndata)
       END IF

       spline = .false.
       wk = 0
       ierr = 0
       iwk = 2*ndata

       CALL PCHEZ(ndata+i_full, xfull, yfull, dyfull, spline,
     1   wk, iwk, ierr)
        IF(ierr.lt.0) STOP 'LEGENDRE: error in SPLINE'

       CALL PCHEV(ndata+i_full, xfull, yfull, dyfull,
     1   npts, x, y, dy, ierr)
        IF(ierr.lt.0) STOP 'LEGENDRE: error in EVAL_SPLINE'

       DEALLOCATE(xfull, yfull, dyfull, wk)

       END SUBROUTINE spline_it
