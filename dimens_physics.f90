    module dimens_physics

    doubleprecision,allocatable :: cod_x(:,:),cod_y(:,:)
    doubleprecision,allocatable :: &        
        gxyf0(:,:),    & 
		uxyf0(:,:),    & 
		sxyf(:,:),    &        
        fxyf(:,:),    & 
		hxyf(:,:),    & 
		gxyf(:,:),    &
        uxyf(:,:)
        
    
    doublecomplex,allocatable ::    &
        hxyf_DFT(:,:),    &
        fxyf_DFT(:,:),    &
        gxyf_DFT(:,:),    &
        uxyf_DFT(:,:)
    end module