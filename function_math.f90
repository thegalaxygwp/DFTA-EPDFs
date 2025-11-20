    module function_math
    use allmodules
    implicit doubleprecision (a-h,o-z)
    contains
    !=====================================================================!   
    subroutine function_initial(i,j,temp_x,temp_y,gsf_point)
    integer :: m,i,j
    doubleprecision :: gsf_point(9)
    doubleprecision qxf,qxf_x,qxf_xx
    doubleprecision pyf,pyf_y,pyf_yy,ssxyf
    doubleprecision sgxyf,sgxyf_x,sgxyf_xx,sgxyf_y,sgxyf_yy,sfxyf
    doubleprecision suxyf,suxyf_x,suxyf_xx,suxyf_y,suxyf_yy,shxyf

    temp_L  =(temp_x-half_L)/exp_L1      
	
    qxf     =   1*sin(2*pi*temp_x)*exp(-temp_L**2)
    qxf_x   =   2*exp(-temp_L**2)*(pi*cos(2*pi*temp_x)-temp_L/exp_L1*sin(2*pi*temp_x))
    qxf_xx  =   2*exp(-temp_L**2)*((2*temp_L**2/exp_L1**2-1/exp_L1**2-2*pi**2)*sin(2*pi*temp_x)-4*pi*temp_L/exp_L1*cos(2*pi*temp_x))
    
    !rxf     =   1*cos(2*pi*temp_x)*exp(-temp_L**2)
    !rxf_x   =   2*exp(-temp_L**2)*(-pi*sin(2*pi*temp_x)-temp_L/exp_L1*cos(2*pi*temp_x))
    !rxf_xx  =   2*exp(-temp_L**2)*((2*temp_L**2/exp_L1**2-1/exp_L1**2-2*pi**2)*cos(2*pi*temp_x)+4*pi*temp_L/exp_L1*sin(2*pi*temp_x))
    
	!rxf     =   cos(2d0*pi*temp_x)*tanh(temp_L)
 !   rxf_x   =   -2d0*pi*SIN(2d0*pi*temp_x)*TANH(temp_L)+cos(2d0*pi*temp_x)*(1-tanh(temp_L)**2)/exp_L1
 !   rxf_xx  =   -4d0*pi**2*cos(2d0*pi*temp_x)*TANH(temp_L)-4d0*pi*sin(2d0*pi*temp_x)*(1-tanh(temp_L)**2)/exp_L1-2d0*cos(2d0*pi*temp_x)*tanh(temp_L)*(1-tanh(temp_L)**2)/exp_L1**2d0
     
	temp_L  =   (temp_x-half_L)/exp_L2
	rxf     =   tanh(temp_L)
    rxf_x   =   (1-tanh(temp_L)**2)/exp_L2
    rxf_xx  =   -2d0*tanh(temp_L)*(1-tanh(temp_L)**2)/exp_L2**2d0
	
	if(temp_x <=0.1)then	
		temp_x =0.1
		temp_L  =(temp_x-half_L)/exp_L2      
		rxf     =   tanh(temp_L)
	    rxf_x   =0
	    rxf_xx  =0
	elseif(temp_xi>=0.9)then
		temp_x =0.9
		temp_L  =(temp_x-half_L)/exp_L2      
		rxf     =   tanh(temp_L)
	    rxf_x   =0
	    rxf_xx  =0
		
	endif

	temp_L  =   (temp_x-half_L)/exp_L3
    !ssxyf     =  1d0*(1d0+0.5d0*exp(-temp_L**2)*exp(-((temp_y-half_L)/exp_L1)**2))!*pyf
    ssxyf     =  1d0+0.1d0*tanh(temp_L)
	!ssxyf	  =	 -ssxyf

    pyf     =0;
    pyf_y   =0;
    pyf_yy  =0;    
    do m=1,1
        pyf     =   pyf+cos(2*pi*m*temp_y)
        pyf_y   =   pyf_y-2*m*pi*sin(2*m*pi*temp_y)        
        pyf_yy  =   pyf_yy-(2*m*pi)**2*cos(2*m*pi*temp_y)
	enddo
    
	!temp_L  =(temp_y-half_L)/exp_L1      
	!
 !   pyf     =   1*sin(2*pi*temp_y)*exp(-temp_L**2)
 !   pyf_y   =   2*exp(-temp_L**2)*(pi*cos(2*pi*temp_y)-temp_L/exp_L1*sin(2*pi*temp_y))
 !   pyf_yy  =   2*exp(-temp_L**2)*((2*temp_L**2/exp_L1**2-1/exp_L1**2-2*pi**2)*sin(2*pi*temp_y)-4*pi*temp_L/exp_L1*cos(2*pi*temp_y))
	
    sgxyf     =  qxf*pyf
    sgxyf_x   =  qxf_x*pyf
    sgxyf_xx  =  qxf_xx*pyf
    sgxyf_y   =  qxf*pyf_y
    sgxyf_yy  =  qxf*pyf_yy
		
    suxyf     =  rxf*pyf
    suxyf_x   =  rxf_x*pyf
    suxyf_xx  =  rxf_xx*pyf
    suxyf_y   =  rxf*pyf_y
    suxyf_yy  =	 rxf*pyf_yy
    
    
    sfxyf     =  sgxyf-ssxyf*(sgxyf_xx+sgxyf_yy)
    shxyf	  =  suxyf-ssxyf*(suxyf_xx+suxyf_yy)
    
	gsf_point(1)    = sgxyf
	gsf_point(2)    = suxyf
	gsf_point(3)    = sfxyf
	gsf_point(4)    = shxyf
	gsf_point(5)    = ssxyf	
	
    end subroutine function_initial

    end module function_math