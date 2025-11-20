    subroutine transprt()
    use allmodules
    use extend_data
    implicit doubleprecision (a-h,o-z)
    !--------simulating parameters-----------------!
    nreg     =1
    nrecord  =0
    !---------math parameters----------------------!
    fw_part     =(1q0,0q0)
    bh_part     =(0q0,1q0)
    pi          =3.1415926535897932384626433832795028841971693993751q0
    de          =1.d0
    sqrt_2      =sqrt(2q0)
    kappa       =1d0        !'n_ka==0 , for uniform sxyf condition'
    if(kappa==0)then
        n_ka    =0
    else
        n_ka    =1
    endif    
    n_reduce    =0             !    use 0 while for linear sigle mode
    n_dimens    =2
    n_ctp       =1
    solution_n  =1
    mode_N      =1			   !    the sigle num of y mode for linear type
    if(n_reduce==0) then
		nlr     =0
	elseif(n_reduce/=0) then
		nlr     =1
	endif	
    !---------coordinate parameters----------------!
    nj          =32           !    if(nlr==1) nj must be bigger e.g. 64
    r_a         =1.d0          !    width of x
    exp_L1      =0.1d0         !    function parameter	
    exp_L2      =0.02d0        !    function parameter	
    exp_L3      =0.2d0        !    function parameter
    half_L      =0.5d0*r_a     !    function parameter 
    r_in        =0d0           !    inner boundary of x
    r_ot        =r_in+r_a      !    inner boundary of x
    r_t         =1.d0          !    *r_a
    dy          =r_t/nj        !    *dy
    nd          =5          !    grid order of x ,>>ni
    ni          =2**nd+1       !    ni
    dx          =r_a/(ni-1)    !    dx
    dg_s        =r_t
    dg_r        =2q0*r_a
    dg_mode     =r_t/real(mode_N)
	nidi        =5
	njdj        =nj/4
    !---------setting coordinate array-------------!
    nr          =2**(nd+1)     !    >>ni
    ns          =nj            !    >>nj
    nbi         =(ni-2)*nj     !    Dimens Of A Matrix Without Boundary
    nbj         =nbi           !    Dimens Of A Matrix

    nyi         =ni*nj         !    Dimens Of A Matrix With Boundary
    nyj         =nyi           !    Dimens Of A Matrix
    kbi         =(ni-2)*nj     !    Dimens Of A Matrix
    kbj         =kbi           !    Dimens Of A Matrix

    lnyi        =ni*2          !    Dimens Of Fast Solver
    lnyj        =lnyi          !    Dimens Of Fast Solver
    lnbi        =(ni-2)*2      !    Dimens Of Fast Solver
    lnbj        =lnbi          !    Dimens Of Fast Solver
    !---------setting parallel thread -------------!
    max_thread  =32          !    Manually Set The Number Of Parallel Threads 
    write(*,*)'we use',max_thread,'threads'
	write(*,*)'x_grid',ni
	write(*,*)'y_grid',nj

    Nj_th   =max_thread
    if(ns<max_thread)  Nj_th   =ns

    return
    end