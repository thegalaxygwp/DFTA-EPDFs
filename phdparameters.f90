module phdparameters
	implicit doubleprecision (a-h,o-z)
    include 'omp_lib.h'
    
    
	doubleprecision :: pi,sqrt_2 !math parameters!
	doublecomplex :: fw_part,bh_part  !math parameters!
	integer :: nd,ni,nj,nf,nbi,nbj,nyi,nyj,solution_n,mode_N,case_N,nidi,njdj  !number of field!
    integer :: lnyi,lnyj,lnbi,lnbj
	integer :: nstep,nt,irecord,nwrite,nreg,nrecord
    integer :: nr,ns,nlr,linear_case,kbi,kbj
    integer :: n_ka,n_reduce,n_dimens,n_ctp,nc_stop,ncut_en
    integer :: num_core,num_thread,max_thread,Nr1_dy,Nr2_dy,N1_dy,N2_dy,Nj_th,Ns_dy
    Integer , parameter :: FFT_Forward = -1
    Integer , parameter :: FFT_Inverse = 1
    integer :: nport1,nport2,nport3,nport4,nport5,nport6
    doubleprecision :: dg_r,dg_s,dg_mode
	doubleprecision :: cfl,t,dt,rate_t,rate_b,rate_k,rate_f,dt_0      !simulating parameters!
	doubleprecision :: epsilon,eta,alpha,de,dethe,dp_dx !phusical parameters!
	doubleprecision :: by_0,bz_0,bt_0,pe_0,re_0,beta0,kappa,gamma !normalized parameters!
	doubleprecision :: r_a,r_t,r_in,r_ot,half_L,dx,dy,dr1,ds1,dr2,ds2  !coordinate parameters!
	doubleprecision :: amp,exp_L1,a2_loc,a_amp,exp_L2,exp_L3  !perturbation parameters!
	doubleprecision :: energy_k,energy_bx,energy_p
    doubleprecision :: h1,h2,h3,h4
    doubleprecision :: w_drift,o_point,o_point1,o_point2,x_point
    doubleprecision :: bxt_q,bxt_h !for analysis
    real :: tstart,tend
    
    character*128 filedata
    doubleprecision :: ffmax1,ffmax2,ffmax3,ffmax4,ffmax5,ffmax6,ffmax7
    doubleprecision :: carryp,fixbz_0,beta_c

    
	end module