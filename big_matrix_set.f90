    subroutine big_matrix_set()
    use allmodules
    implicit doubleprecision (a-h,o-z)
    
    nkr=2
    if(n_ka/=0)then     !for non-uniform sxyf condition
        print*, 'n_ka/=0 , use Algbra FFT algorithm'
        if(n_reduce==0)then ! m_mode limited
            print*, 'use Fast Matching Algbra FFT algorithm'
            allocate(Temp_rsin(lnbi),Temp_Xxz(lnbi),linear_sps(nkr))
            allocate(Bmgt_ma1(nr,ns),Bmgt_ma2(nr,ns))
            linear_sps=(/mode_N+1,ns-mode_N+1/)
            call big_matrix_x_fast()
            call big_matrix_y_fast()
        elseif(n_reduce/=0)then  ! m_mode nonlinear
            print*, 'use normal Algbra FFT algorithm'
            allocate(Temp_sin(nbi),Temp_Xxz(nbi),linear_sps(nkr))
            linear_sps=(/mode_N+1,ns-mode_N+1/)
            call big_matrix_x()
            call big_matrix_y()
        endif
    else
        print*, 'n_ka==0 , for uniform sxyf condition'        
    endif
        
    return
    end