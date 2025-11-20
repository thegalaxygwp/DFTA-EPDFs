    module dimens_solver
    doublecomplex,allocatable :: X_T(:,:),XL_T1(:,:),XL_T2(:,:)
    doublecomplex,allocatable :: Y_T(:,:),YL_T1(:,:),YL_T2(:,:)
    doublecomplex,allocatable :: Z_T(:,:),ZL_T1(:,:),ZL_T2(:,:)
    
    doublecomplex,allocatable :: Bmgt_ma1(:,:),Bmgt_ma2(:,:)
    
    doublecomplex,allocatable :: Temp_cos(:),Temp_Xyy(:)
    doublecomplex,allocatable :: Temp_sin(:),Temp_Xxz(:)
    
    doublecomplex,allocatable :: Temp_rcos(:),Temp_rXyy(:)
    doublecomplex,allocatable :: Temp_rsin(:),Temp_rXxz(:)
        
    integer :: nkr
    integer,allocatable :: linear_sps(:)
    
    end module