module extend_data
	implicit none
	contains
	
	real(8)function real_4ex8(x)
	real(4) :: x
	real(8) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	real_4ex8=x1
	end function real_4ex8
		
	real(16)function real_4ex16(x)
	real(4)  :: x
	real(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	real_4ex16=x1
	end function real_4ex16
	
	real(16)function real_8ex16(x)
	real(8) :: x
	real(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	real_8ex16=x1
	end function real_8ex16
	
	complex(8)function cmplx_4ex8(x)
	complex(4) :: x
	complex(8) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	cmplx_4ex8=x1
	end function cmplx_4ex8
	
	complex(16)function cmplx_4ex16(x)
	complex(4) :: x
	complex(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	cmplx_4ex16=x1
	end function cmplx_4ex16
	
	complex(16)function cmplx_8ex16(x)
	complex(8) :: x
	complex(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	cmplx_8ex16=x1
    end function cmplx_8ex16
    
    
	complex(16)function cmplx_all_ex16(x)
    complex,intent(in) :: x
	complex(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	cmplx_all_ex16=x1
    end function cmplx_all_ex16
	
    
	real(16)function real_all_ex16(x)
	real,intent(in) :: x
	real(16) :: x1	
	character*1024 exc
	write(exc,*)x
	read(exc,*)x1
	real_all_ex16=x1
    end function real_all_ex16
    
	endmodule extend_data