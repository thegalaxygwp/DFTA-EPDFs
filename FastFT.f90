    Module FastFT
    use allmodules
    Implicit doubleprecision(a-h,o-z)
    contains

    Subroutine fcFFT( x , forback )    !//Subroutine FFT , Cooley-Tukey , radix-2    !// www.fcode.cn
    doublecomplex :: x(:),tempx
    Integer , Intent(IN) :: forback
    Integer :: n
    integer :: i , j , k , ncur , ntmp , itmp
    doubleprecision :: dpi,f_n
    doublecomplex :: ctmp,epi
    n = size(x)
    ncur = n
    dpi=3.1415926535897932384626433832795028841971693993751q0
    Do
        ntmp = ncur
        !epi = 2.0d0 * pi / ncur
        f_n  =1.d0/ncur
        epi  = cmplx(0.d0,2.d0*dpi*f_n)
        ncur = ncur / 2
        if ( ncur < 1 ) exit
        Do j = 1 , ncur
            Do i = j , n , ntmp
                itmp = i + ncur
                ctmp = x(i) - x(itmp)
                x(i) = x(i) + x(itmp)
                !x(itmp) = ctmp * exp( forback * cmplx( 0.0d0 , epi*(j-1) ) )
                x(itmp) = ctmp * exp( forback * epi * (j-1) )
            End Do
        End Do
    End Do
    j = 1
    Do i = 1, n - 1
        If ( i < j ) then
            ctmp = x(j)
            x(j) = x(i)
            x(i) = ctmp
        End If
        k = n/2
        Do while( k < j )
            j = j - k
            k = k / 2
        End Do
        j = j + k
    End Do
    Return
    End Subroutine fcFFT
    !----------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    Subroutine fcFFT2( xx, nnr, nns, forback )
    integer nnr,nns
    doublecomplex :: xx(nnr,nns)
    doublecomplex,allocatable ::x(:),y(:)
    Integer , Intent(IN) :: forback
    integer :: i , j
    allocate(x(nnr))
    allocate(y(nns))
    !$omp parallel shared(xx) num_threads(max_thread)
    !$omp do collapse(1) private(i,y)! schedule(dynamic,Nr2_dy)
    do i=1,nnr
        y(:)=xx(i,:)
        call fcFFT( y , forback )
        xx(i,:)=y(:)
    enddo
    !$omp enddo
    !$omp end parallel
    !$omp parallel shared(xx) num_threads(Nj_th)
    !$omp do collapse(1) private(j,x)!  schedule(dynamic,Ns_dy)    
    do j=1,nns
        x(:)=xx(:,j)
        call fcFFT( x , forback )
        xx(:,j)=x(:)
    enddo
    !$omp enddo
    !$omp end parallel    
    if(forback==FFT_Inverse)then
        !$omp parallel shared(xx)  num_threads(max_thread)
        !$omp do collapse(2) private(i,j)
        do j=1,nns
            do i=1,nnr
                xx(i,j)=xx(i,j)/nnr/nns
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endif
    deallocate(x)
    deallocate(y)
    Return
    End Subroutine fcFFT2
    !----------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    Subroutine dcFFT( x , forback )
    complex(16) :: x(:)
    Integer , Intent(IN) :: forback
    Integer :: n
    integer :: i , j , k , ncur , ntmp , itmp
    real(16) :: dpi,f_n
    complex(16) :: ctmp,epi
    n = size(x)
    ncur = n
    dpi=3.1415926535897932384626433832795028841971693993751q0
    Do
        ntmp = ncur
        !epi = 2.0d0 * pi / ncur
        f_n  =1.q0/ncur
        epi  = qcmplx(0.q0,2.q0*dpi*f_n)
        ncur = ncur / 2
        if ( ncur < 1 ) exit
        Do j = 1 , ncur
            Do i = j , n , ntmp
                itmp = i + ncur
                ctmp = x(i) - x(itmp)
                x(i) = x(i) + x(itmp)
                !x(itmp) = ctmp * exp( forback * cmplx( 0.0d0 , epi*(j-1) ) )
                x(itmp) = ctmp * exp( forback * epi * (j-1) )
            End Do
        End Do
    End Do
    j = 1
    Do i = 1, n - 1
        If ( i < j ) then
            ctmp = x(j)
            x(j) = x(i)
            x(i) = ctmp
        End If
        k = n/2
        Do while( k < j )
            j = j - k
            k = k / 2
        End Do
        j = j + k
    End Do
    Return
    End Subroutine dcFFT
    !----------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    Subroutine dcFFT2( xx, nnr, nns, forback )
    integer nnr,nns
    complex(16) :: xx(:,:)
    complex(16),allocatable ::x(:),y(:)
    Integer , Intent(IN) :: forback
    integer :: i , j
    allocate(x(nnr))
    allocate(y(nns))
    !$omp parallel shared(xx) num_threads(max_thread)
    !$omp do collapse(1) private(i,y)! schedule(dynamic,Nr2_dy)
    do i=1,nnr
        y(:)=xx(i,:)
        call dcFFT( y , forback )
        xx(i,:)=y(:)
    enddo
    !$omp enddo
    !$omp end parallel
    !$omp parallel shared(xx) num_threads(Nj_th)
    !$omp do collapse(1) private(j,x)!  schedule(dynamic,Ns_dy)
    do j=1,nns
        x(:)=xx(:,j)
        call dcFFT( x , forback )
        xx(:,j)=x(:)
    enddo
    !$omp enddo
    !$omp end parallel    
    if(forback==FFT_Inverse)then
        !$omp parallel shared(xx)  num_threads(max_thread)
        !$omp do collapse(2) private(i,j)
        do j=1,nns
            do i=1,nnr
                xx(i,j)=xx(i,j)/nnr/nns
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endif
    deallocate(x)
    deallocate(y)
    Return
    End Subroutine dcFFT2
    !----------------------------------------------------------------------------!

    complex(16) function tranf1(i,dg_d,x)
    integer :: i
    real(16) :: dg_d
    complex(16) :: x
    tranf1=qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*x
    end function tranf1
    
    complex(16) function tranf2(i,dg_d,x)
    integer :: i
    real(16) :: dg_d
    complex(16) :: x
    tranf2=qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*x
    end function tranf2

    complex(16) function tranf3(i,dg_d,x)
    integer :: i
    real(16) :: dg_d
    complex(16) :: x
    tranf3=qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*qcmplx(0.q0,(FFT_Inverse*2.q0*pi*i/dg_d))*x
    end function tranf3
    
    
    !--------------------------------------------------------------

    doubleprecision function sum_inverse(temp_y,ffcns)
    doubleprecision :: temp_y
    doublecomplex :: temp_ans
    doublecomplex :: ffcns(ns)
    temp_ans=(0d0,0d0)
    do i=1,ns
        if(i<=ns/2+1)then
            temp_ans=temp_ans+exp(cmplx(0d0,FFT_Inverse*2d0*pi*temp_y*(i-1)/dg_s))*ffcns(i)/ns
        else
            temp_ans=temp_ans+exp(cmplx(0d0,-FFT_Inverse*2d0*pi*temp_y*(ns+1-i)/dg_s))*ffcns(i)/ns
        endif
    enddo
    sum_inverse=real(temp_ans)
    end function sum_inverse
    
    End Module FastFT
