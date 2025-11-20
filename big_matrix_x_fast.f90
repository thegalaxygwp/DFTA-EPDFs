    subroutine big_matrix_x_fast()
    use cmi_module
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,jk,k,m,n,mi,mj,mk,info
    integer,allocatable :: ipiv1(:)
    doublecomplex,allocatable :: work1(:)
    !=============================================================================================!
    allocate(XL_T1(lnbi,lnbj),XL_T2(lnbi,lnbj))
    !=============================================================================================!
    allocate(ipiv1(lnbi),work1(lnbj))
	!Averaging after two solutions using staggered grids in linear mode
    !=============================================================================================!
    m=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk
            do while(j>ns)
                j=j-ns
            enddo
            n=1
            do mi=2,ni-1
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do k=1,nkr
                    mj=linear_sps(k)
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    ep6= 1d0-dd_t2(mi,mj)*sxyf(i,j)
                    ep7= 1d0-dd_t2(mk,mj)*sxyf(i,j)
                    XL_T1(m,n)=(ep1*ep6-ep2*ep7)*ep3
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo
    !=============================================================================================!
    m=1
    do i=2,ni-1
        do jk=1,nkr
            j=nidi*i+njdj*jk+ns/2
            do while(j>ns)
                j=j-ns
            enddo
            n=1
            do mi=2,ni-1
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do k=1,nkr
                    mj=linear_sps(k)
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    ep6= 1d0-dd_t2(mi,mj)*sxyf(i,j)
                    ep7= 1d0-dd_t2(mk,mj)*sxyf(i,j)
                    XL_T2(m,n)=(ep1*ep6-ep2*ep7)*ep3
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo

    write(*,*)'loading matrix for fast_solve_x system...'
    
    !call zgetrf( lnbi, lnbj, XL_T1, lnbi, ipiv1, info )
    !write(*,*)info
    call zgetrf_hand( lnbi, lnbj, XL_T1, lnbi, ipiv1, info )
    write(*,*)info
    !call zgetri( lnbi, XL_T1, lnbi, ipiv1, work1, lnbj, info )
    !write(*,*)info
    call zgetri_hand( lnbi, XL_T1, lnbi, ipiv1, work1, lnbj, info )
    write(*,*)info
    
    
    !call zgetrf( lnbi, lnbj, XL_T2, lnbi, ipiv1, info )
    !write(*,*)info
    call zgetrf_hand( lnbi, lnbj, XL_T2, lnbi, ipiv1, info )
    write(*,*)info
    !call zgetri( lnbi, XL_T2, lnbi, ipiv1, work1, lnbj, info )
    !write(*,*)info
    call zgetri_hand( lnbi, XL_T2, lnbi, ipiv1, work1, lnbj, info )
    write(*,*)info
    
    deallocate(ipiv1,work1)
    
    write(*,*)'Done!'
    
    return
    end