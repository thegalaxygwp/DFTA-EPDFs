    subroutine big_matrix_x()
    use cmi_module
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,nk,m,n,mi,mj,mk,info
    integer,allocatable :: ipiv1(:)
    doublecomplex,allocatable :: work1(:)
	integer(4) :: status
    !=============================================================================================!
    allocate(X_T(nbi,nbj))
    !=============================================================================================!
    allocate(ipiv1(nbi),work1(nbj))
    !=============================================================================================!
    m=1
    do i=2,ni-1
        do j=1,nj
            n=1
            do mi=2,ni-1
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do mj=1,ns
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    ep4= 1d0-dd_t2(mi,mj)*sxyf(i,j)
                    ep5= 1d0-dd_t2(mk,mj)*sxyf(i,j)
                    X_T(m,n)=(ep1*ep4-ep2*ep5)*ep3
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo
    
    write(*,*)'loading matrix for x system...'
	
    call mkl_set_num_threads(max_thread)
    !call zgetrf( nbi, nbj, X_T, nbi, ipiv1, info )
    call zgetrf_hand( nbi, nbj, X_T, nbi, ipiv1, info )
    write(*,*)info
    !
    !call zgetri( nbi, X_T, nbi, ipiv1, work1, nbj, info )
	call zgetri_hand( nbi, X_T, nbi, ipiv1, work1, nbj, info )
    write(*,*)info
    
    !call complex_matrix_inverse_omp_di(X_T, nbi)
    
    !call complex_matrix_inverse_omp(X_T, nbi)
    
    deallocate(ipiv1,work1)
    write(*,*)'Done!'
    
    return
    end