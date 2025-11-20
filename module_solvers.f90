module module_solvers
	
    use FastFT
    contains    
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_x_fast_linear(fftxxx,fddxxx)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j,k,jk,ik   
    doubleprecision ::  fddxxx(:,:)
    doublecomplex ::  fftxxx(:,:)
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rsin(k)=nr*ns*cmplx(fddxxx(i,j),0d0)
            k=k+1
        enddo
    enddo

    Temp_rXxz=matmul(XL_T1,Temp_rsin)

    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma1(i,j)=Temp_rXxz(k)
            k=k+1
        enddo
    enddo
    Bmgt_ma1(1,:)=(0d0,0d0)
    Bmgt_ma1(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma1(i,j)=-Bmgt_ma1(nr+2-i,j)
        enddo
    enddo
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk+ns/2
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rsin(k)=nr*ns*cmplx(fddxxx(i,j),0d0)
            k=k+1
        enddo
    enddo

    Temp_rXxz=matmul(XL_T2,Temp_rsin)

    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma2(i,j)=Temp_rXxz(k)
            k=k+1
        enddo
    enddo
    Bmgt_ma2(1,:)=(0d0,0d0)
    Bmgt_ma2(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma2(i,j)=-Bmgt_ma2(nr+2-i,j)
        enddo
    enddo
    !=============================================================================================!
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            fftxxx(i,j)=0.5d0*(Bmgt_ma1(i,j)+Bmgt_ma2(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel    
    return
    end subroutine solve_x_fast_linear
    
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!    
    subroutine solve_y_fast_linear(fftyyy,fddyyy)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j,k,jk,ik    
    doubleprecision  fddyyy(:,:)
    doublecomplex  fftyyy(:,:)
    !=============================================================================================!
    k=1    
    do i=1,ni
        do jk=1,nkr
            j=nidi*i+njdj*jk
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rcos(k)=nr*ns*cmplx(fddyyy(i,j),0d0)
            k=k+1
        enddo
    enddo
    
    Temp_rXyy=matmul(YL_T1,Temp_rcos)
    
    k=1
    do i=1,ni
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma1(i,j)=Temp_rXyy(k)
            k=k+1
        enddo
    enddo    
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma1(i,j)=Bmgt_ma1(nr+2-i,j)
        enddo
    enddo
    !-----------------------------------------------------------------------------------!
    k=1    
    do i=1,ni
        do jk=1,nkr
            j=nidi*i+njdj*jk+ns/2
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rcos(k)=nr*ns*cmplx(fddyyy(i,j),0d0)
            k=k+1
        enddo
    enddo
    
    Temp_rXyy=matmul(YL_T2,Temp_rcos)
    
    k=1
    do i=1,ni
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma2(i,j)=Temp_rXyy(k)
            k=k+1
        enddo
    enddo
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma2(i,j)=Bmgt_ma2(nr+2-i,j)
        enddo
    enddo
    !=============================================================================================!
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            fftyyy(i,j)=0.5d0*(Bmgt_ma1(i,j)+Bmgt_ma2(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
    end subroutine solve_y_fast_linear
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_z_fast_linear(fftzzz,fddzzz)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j,k,jk,ik    
    doubleprecision  fddzzz(:,:)
    doublecomplex  fftzzz(:,:)
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rsin(k)=nr*ns*cmplx(fddzzz(i,j),0d0)
            k=k+1
        enddo
    enddo

    Temp_rXxz=matmul(ZL_T1,Temp_rsin)

    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma1(i,j)=Temp_rXxz(k)
            k=k+1
        enddo
    enddo
    Bmgt_ma1(1,:)=(0d0,0d0)
    Bmgt_ma1(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma1(i,j)=-Bmgt_ma1(nr+2-i,j)
        enddo
    enddo
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk+ns/2
            do while(j>ns)
                j=j-ns
            enddo
            Temp_rsin(k)=nr*ns*cmplx(fddzzz(i,j),0d0)
            k=k+1
        enddo
    enddo

    Temp_rXxz=matmul(ZL_T2,Temp_rsin)

    k=1
    do i=2,ni-1
        do jk=1,nkr
            j=linear_sps(jk)
            Bmgt_ma2(i,j)=Temp_rXxz(k)
            k=k+1
        enddo
    enddo
    Bmgt_ma2(1,:)=(0d0,0d0)
    Bmgt_ma2(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            Bmgt_ma2(i,j)=-Bmgt_ma2(nr+2-i,j)
        enddo
    enddo
    !=============================================================================================!
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            fftzzz(i,j)=0.5d0*(Bmgt_ma1(i,j)+Bmgt_ma2(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel    
    return
    end subroutine solve_z_fast_linear
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_x(fftxxx,fddxxx)
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,m,n,mi,mj,mk,info
    doubleprecision  fddxxx(:,:)
    doublecomplex  fftxxx(:,:)
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do j=1,nj
            Temp_sin(k)=nr*ns*cmplx(fddxxx(i,j),0d0)
            k=k+1
        enddo
    enddo
    Temp_Xxz=matmul(X_T,Temp_sin)
    do k=1,nbi
        i=k/nj+2
        j=mod(k,nj)
        if(j==0)j=nj
        if(j==nj)i=i-1
        fftxxx(i,j)=Temp_Xxz(k)
    enddo
    fftxxx(1,:)=(0d0,0d0)
    fftxxx(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            fftxxx(i,j)=-fftxxx(nr+2-i,j)
        enddo
    enddo    
    return
    end subroutine solve_x
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_y(fftyyy,fddyyy)
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,m,n,mi,mj,mk,info
    doubleprecision  fddyyy(:,:)
    doublecomplex  fftyyy(:,:)
    !=============================================================================================!
    k=1
    do i=1,ni
        do j=1,nj
            Temp_cos(k)=nr*ns*cmplx(fddyyy(i,j),0d0)
            k=k+1
        enddo
    enddo
    Temp_Xyy=matmul(Y_T,Temp_cos)
    do k=1,nyi
        i=k/nj+1
        j=mod(k,nj)
        if(j==0)j=nj
        if(j==nj)i=i-1
        fftyyy(i,j)=Temp_Xyy(k)
    enddo
    do i=ni+1,nr
        do j=1,ns
            fftyyy(i,j)=fftyyy(nr+2-i,j)
        enddo
    enddo
    return
    end subroutine solve_y
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_z(fftzzz,fddzzz)
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,m,n,mi,mj,mk,info
    doubleprecision  fddzzz(:,:)
    doublecomplex  fftzzz(:,:)
    !=============================================================================================!
    k=1
    do i=2,ni-1
        do j=1,nj
            Temp_sin(k)=nr*ns*cmplx(fddzzz(i,j),0d0)
            k=k+1
        enddo
    enddo
    Temp_Xxz=matmul(Z_T,Temp_sin)
    do k=1,nbi
        i=k/nj+2
        j=mod(k,nj)
        if(j==0)j=nj
        if(j==nj)i=i-1
        fftzzz(i,j)=Temp_Xxz(k)
    enddo
    fftzzz(1,:)=(0d0,0d0)
    fftzzz(ni,:)=(0d0,0d0)
    do i=ni+1,nr
        do j=1,ns
            fftzzz(i,j)=-fftzzz(nr+2-i,j)
        enddo
    enddo
    return
    end subroutine solve_z
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_uniform_cos(f00cos,fftcos)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doublecomplex :: fftcos(:,:),f00cos(:,:)
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            f00cos(i,j)=fftcos(i,j)/de_t2(i,j)
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine solve_uniform_cos
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine solve_uniform_sin(f00sin,fftsin)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doublecomplex :: fftsin(:,:),f00sin(:,:)
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            f00sin(i,j)=fftsin(i,j)/de_t2(i,j)
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine solve_uniform_sin
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine trans_sin(fftsin,fddsin)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doubleprecision :: fddsin(:,:)
    doublecomplex :: fftsin(:,:)
    !=============================================================================================!    
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            if(i<=ni)then
                fftsin(i,j)= cmplx(fddsin(i,j),0d0)
            else
                fftsin(i,j)=-cmplx(fddsin(nr+2-i,j),0d0)
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    call fcFFT2( fftsin, nr, ns, FFT_Forward )
    fftsin(1,:) =cmplx(0d0,0d0)
    fftsin(ni,:)=cmplx(0d0,0d0)
    call fftsfelter(fftsin)
    return
    end subroutine trans_sin
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine trans_cos(fftcos,fddcos)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doubleprecision :: fddcos(:,:)
    doublecomplex :: fftcos(:,:)
    !=============================================================================================!    
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,ns
        do i=1,nr
            if(i<=ni)then
                fftcos(i,j)= cmplx(fddcos(i,j),0d0)
            else
                fftcos(i,j)= cmplx(fddcos(nr+2-i,j),0d0)
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    call fcFFT2( fftcos, nr, ns, FFT_Forward )
    call fftsfelter(fftcos)
    return
    end subroutine trans_cos
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine inves_sin(fftsin,fddsin)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doubleprecision :: fddsin(:,:)
    doublecomplex :: fftsin(:,:)
    doublecomplex,allocatable :: fiisin(:,:)
    allocate(fiisin(nr,ns))
    fiisin=fftsin
    call fcFFT2( fiisin, nr, ns, FFT_Inverse )
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,nj
        do i=1,ni
            fddsin(i,j)= real(fiisin(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(fiisin)
    fddsin(1,:) =0d0
    fddsin(ni,:)=0d0
    return
    end subroutine inves_sin
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!
    subroutine inves_cos(fftcos,fddcos)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doubleprecision :: fddcos(:,:)
    doublecomplex :: fftcos(:,:)
    doublecomplex,allocatable :: fiicos(:,:)
    allocate(fiicos(nr,ns))
    fiicos=fftcos
    call fcFFT2( fiicos, nr, ns, FFT_Inverse )
    !$omp parallel num_threads(max_thread)
    !$omp do collapse(2) private(i,j)
    do j=1,nj
        do i=1,ni
            fddcos(i,j)= real(fiicos(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(fiicos)
    return
    end subroutine inves_cos
    !=============================================================================================!
    !=============================================================================================!
    !=============================================================================================!    
    subroutine fftsfelter(fftsers)
    implicit doubleprecision (a-h,o-z)
    integer :: i,j
    doublecomplex :: fftsers(:,:)
    if(nlr==0)then
        !$omp parallel num_threads(max_thread)
        !$omp do collapse(2) private(i,j)
        do i=1,nr
            do j=1,ns
                if(j/=2+mode_N-1 .and. j/=ns-mode_N+1 .or. (i>ni-nr/4 .and. i<ni+nr/4))then
                    fftsers(i,j)=cmplx(0d0,0d0)
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    elseif(nlr==1)then        
        !$omp parallel num_threads(max_thread)
        !$omp do collapse(2) private(i,j)
        do i=1,nr
            do j=1,ns
                if((j>ns/2+1-ns/4 .and. j<ns/2+1+ns/4) .or. (i>ni-nr/4 .and. i<ni+nr/4))then
                    fftsers(i,j)=cmplx(0d0,0d0)
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endif
    
    return
    end subroutine fftsfelter
    
    end module module_solvers