    subroutine big_matrix_y_fast()
    use cmi_module
    implicit doubleprecision (a-h,o-z)
    doublecomplex ::ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8
    integer :: i,j,jk,k,m,n,mi,mj,mk,info
    integer,allocatable :: ipiv(:)
    doublecomplex,allocatable :: work(:)
    !=============================================================================================!
    allocate(YL_T1(lnyi,lnyj),YL_T2(lnyi,lnyj),Temp_rcos(lnyi),Temp_rXyy(lnyi))
    !=============================================================================================!
    allocate(ipiv(lnyi),work(lnyj))
    !=============================================================================================!
    m=1
    do i=1,ni
        do jk=1,nkr
            j=nidi*i+njdj*jk
            do while(j>ns)
                j=j-ns
            enddo
            n=1
            do mi=1,ni
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do k=1,nkr
                    mj=linear_sps(k)
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    if(mi==1)then
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        YL_T1(m,n)=ep1*ep4*ep3
                    elseif(mi==ni)then
                        ep4    = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        YL_T1(m,n)=ep1*ep4*ep3
                    else
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        ep5     = 1d0-dd_t2(mk,mj)*sxyf(i,j)
                        YL_T1(m,n)=(ep1*ep4+ep2*ep5)*ep3
                    endif
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo

    !=============================================================================================!
    m=1
    do i=1,ni
        do jk=1,nkr
            j=nidi*i+njdj*jk+ns/2
            do while(j>ns)
                j=j-ns
            enddo
            n=1
            do mi=1,ni
                mk=nr+2-mi
                ep1= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mi-1)*(i-1)/nr))
                ep2= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mk-1)*(i-1)/nr))
                do k=1,nkr
                    mj=linear_sps(k)
                    ep3= exp(dcmplx(0d0,FFT_Inverse*2d0*pi*(mj-1)*(j-1)/ns))
                    if(mi==1)then
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        YL_T2(m,n)=ep1*ep4*ep3
                    elseif(mi==ni)then
                        ep4    = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        YL_T2(m,n)=ep1*ep4*ep3
                    else
                        ep4     = 1d0-dd_t2(mi,mj)*sxyf(i,j)
                        ep5     = 1d0-dd_t2(mk,mj)*sxyf(i,j)
                        YL_T2(m,n)=(ep1*ep4+ep2*ep5)*ep3
                    endif
                    n=n+1
                enddo
            enddo
            m=m+1
        enddo
    enddo


    write(*,*)'loading matrix for fast_solve_y system...'
    call zgetrf( lnyi, lnyj, YL_T1, lnyi, ipiv, info )
    write(*,*)info
    !call zgetri( lnyi, YL_T1, lnyi, ipiv, work, lnyj, info )
    !write(*,*)info
    call zgetri_hand( lnyi, YL_T1, lnyi, ipiv, work, lnyj, info )
    write(*,*)info
    
    call zgetrf( lnyi, lnyj, YL_T2, lnyi, ipiv, info )
    write(*,*)info
    !call zgetri( lnyi, YL_T2, lnyi, ipiv, work, lnyj, info )
    !write(*,*)info
    call zgetri_hand( lnyi, YL_T2, lnyi, ipiv, work, lnyj, info )
    write(*,*)info
    
    deallocate(ipiv,work)
    write(*,*)'Done!'
    
    return
    end