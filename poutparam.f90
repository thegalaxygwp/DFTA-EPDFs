subroutine poutparam()
    use allmodules
    implicit doubleprecision (a-h,o-z)    
    
	open(153,file='parameters.txt',defaultfile=trim(filedata))	
	write(153,135) ni,nj,nreg                  &
		,nlr,kappa                                 & !simulating parameters!
		,de &!phusical parameters!
        ,n_ka,n_reduce,n_ctp &
		,r_a                    &!normalized parameters!
		,r_in,r_ot,half_L,dx,dy   & !coordinate parameters!
		,pi                            & !math parameters!
		,exp_L1                          !perturbation parameters!	
		135	format(25X,'ni      = ',I6/&
				   25X,'nj      = ',I6/&
				   25X,'nreg    = ',I6/&
				   25X,'nlr     = ',i6/&
				   25X,'kappa   = ',F6.4/&
				   25X,'de      = ',F6.4/&
				   25X,'n_ka    = ',I6/&
				   25X,'nreduce = ',I6/&
				   25X,'n_ctp   = ',I6/&
				   25X,'r_a     = ',F6.4/&	
				   25X,'r_in    = ',F6.4/&		
				   25X,'r_ot    = ',F6.4/&			
				   25X,'half_L     = ',F6.4/&			
				   25X,'dx      = ',F6.4/&			
				   25X,'dy      = ',F6.4/&
				   25X,'pi      = ',F6.4/&
				   25X,'exp_L1  = ',F6.4)
	close(153)
	!stop
	return
	end