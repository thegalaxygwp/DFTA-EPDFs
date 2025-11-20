    subroutine dimens_allocate()
    use allmodules    
    !vtt : vorticity
    !vlt : velocity
    !crt : current
    !mgt : magnetic
    !pse : pressure
    !dst : density
    allocate(cod_x(ni,nj),cod_y(ni,nj))
	
	allocate( &
        dx_t1(nr,ns),dx_t2(nr,ns),dx_t3(nr,ns),   &
        dy_t1(nr,ns),dy_t2(nr,ns),dy_t3(nr,ns),   &
        dxy_1(nr,ns),dyy_x(nr,ns),dxx_y(nr,ns),   &
        dd_t2(nr,ns),de_t2(nr,ns) )                          !dimens_oprator_allocate

    allocate( &
        fxyf(ni,nj),hxyf(ni,nj),sxyf(ni,nj),    &
        gxyf0(ni,nj),uxyf0(ni,nj),    &        
        gxyf(ni,nj),uxyf(ni,nj)  )
        
    
    allocate(    &
        hxyf_DFT(nr,ns),    &
        fxyf_DFT(nr,ns),    &
        gxyf_DFT(nr,ns),    &
		uxyf_DFT(nr,ns)     )
    return
    end