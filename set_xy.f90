subroutine set_xy()
    use allmodules
    implicit doubleprecision (a-h,o-z)
    !-----------------------------coordinlate------------------------!
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,nj
        do i=1,ni
			cod_x(i,j)=real(i-1)*dx+r_in
			cod_y(i,j)=real(j-1)*dy!+r_in
        enddo
    enddo
    !$omp end do
    !$omp end parallel    
    !---------------------------matrix for differentiation-----------!
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            if(i<=nr/2+1)then
                dx_t1(i,j)=cmplx(0d0,(FFT_Inverse*2d0*pi*(i-1)/dg_r))
                dx_t2(i,j)=cmplx(0d0,(FFT_Inverse*2d0*pi*(i-1)/dg_r))**2d0
                dx_t3(i,j)=cmplx(0d0,(FFT_Inverse*2d0*pi*(i-1)/dg_r))**3d0
            else
                dx_t1(i,j)= -cmplx(0d0,(FFT_Inverse*2d0*pi*(nr+1-i)/dg_r))
                dx_t2(i,j)=(-cmplx(0d0,(FFT_Inverse*2d0*pi*(nr+1-i)/dg_r)))**2d0
                dx_t3(i,j)=(-cmplx(0d0,(FFT_Inverse*2d0*pi*(nr+1-i)/dg_r)))**3d0
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel    
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            if(j<=ns/2+1)then
                dy_t1(i,j)= cmplx(0d0,(FFT_Inverse*2d0*pi*(j-1)/dg_s))
                dy_t2(i,j)=(+cmplx(0d0,(FFT_Inverse*2d0*pi*(j-1)/dg_s)))**2d0
                dy_t3(i,j)=(+cmplx(0d0,(FFT_Inverse*2d0*pi*(j-1)/dg_s)))**3d0
            else
                dy_t1(i,j)= -cmplx(0d0,(FFT_Inverse*2d0*pi*(ns+1-j)/dg_s))
                dy_t2(i,j)=(-cmplx(0d0,(FFT_Inverse*2d0*pi*(ns+1-j)/dg_s)))**2d0
                dy_t3(i,j)=(-cmplx(0d0,(FFT_Inverse*2d0*pi*(ns+1-j)/dg_s)))**3d0
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            dd_t2(i,j)=dx_t2(i,j)+dy_t2(i,j)
            dxy_1(i,j)=dy_t1(i,j)*dx_t1(i,j)
            dyy_x(i,j)=dy_t2(i,j)*dx_t1(i,j)
            dxx_y(i,j)=dx_t2(i,j)*dy_t1(i,j)
            de_t2(i,j)=(1d0,0d0)-(dx_t2(i,j)+dy_t2(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
    end