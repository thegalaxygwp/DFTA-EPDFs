    subroutine big_matrix_y()
    use cmi_module
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,nk,m,n,mi,mj,mk,info
    integer,allocatable :: ipiv(:)
    doublecomplex,allocatable :: work(:)
	integer(4) :: status
    !=============================================================================================!
    allocate(Y_T(nyi,nyj),Temp_cos(nyi),Temp_Xyy(nyi))
    !=============================================================================================!
    allocate(ipiv(nyi),work(nyj))
    !=============================================================================================!
    m=1
    do i=1,ni
        do j=1,nj
            n=1
            do mi=1,ni
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do mj=1,ns
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    if(mi==1)then
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        Y_T(m,n)=ep1*ep4*ep3
                    elseif(mi==ni)then                        
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        Y_T(m,n)=ep1*ep4*ep3
                    else
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        ep5     = 1d0-dd_t2(mk,mj)*sxyf(i,j)                     
                        Y_T(m,n)=(ep1*ep4+ep2*ep5)*ep3
                    endif
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo        
	
    write(*,*)'loading matrix for Y system...'
	call mkl_set_num_threads(max_thread)
    !call zgetrf( nyi, nyj, Y_T, nyi, ipiv, info )
    call zgetrf_hand( nyi, nyj, Y_T, nyi, ipiv, info )
    write(*,*)info
	
    !call zgetri( nyi, Y_T, nyi, ipiv, work, nyj, info )
	
    call zgetri_hand( nyi, Y_T, nyi, ipiv, work, nyj, info )
    write(*,*)info
    
    deallocate(ipiv,work)
    write(*,*)'Done!'

    return
    end