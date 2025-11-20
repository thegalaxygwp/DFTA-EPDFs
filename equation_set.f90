    subroutine equation_set()
    use function_math
    implicit doubleprecision (a-h,o-z)
    doubleprecision temp_x,temp_y
    doubleprecision :: gsf_point(9)
    !-----------------------------caculate ------------------------!
    
    do i=1,ni
        do j=1,nj
            temp_x  =cod_x(i,j)
            temp_y  =cod_y(i,j)
			!if(i ==1 )then
			!	temp_x  =cod_x(2,j)
			!elseif(i==ni)then
			!	temp_x  =cod_x(ni-1,j)
			!endif			
            call function_initial(i,j,temp_x,temp_y,gsf_point)
			            
            gxyf0(i,j)= gsf_point(1)
            uxyf0(i,j)= gsf_point(2)	
			fxyf(i,j) = gsf_point(3)
            hxyf(i,j) = gsf_point(4)
			sxyf(i,j) = gsf_point(5)
						
            gxyf(i,j) = 0
            uxyf(i,j) = 0
        enddo
	enddo   
	
	
	open(177,file='function_data0.dat',defaultfile=trim(filedata))
    write(177,*)'x y gxyf0 uxyf0 fxyf hxyf sxyf gxyf uxyf'
    do i=1,ni
        do j=1,nj+1
            if(j<=nj)then
                write(177,1312) cod_x(i,j),cod_y(i,j), &
                    gxyf0(i,j),uxyf0(i,j), &
                    fxyf(i,j) ,hxyf(i,j) , &
					sxyf(i,j) ,gxyf(i,j) , uxyf(i,j)
            elseif(j==nj+1)then
                write(177,1312) cod_x(i,1),cod_y(i,1), &
                    gxyf0(i,1),uxyf0(i,1), &
                    fxyf(i,1) ,hxyf(i,1) , &
					sxyf(i,1) ,gxyf(i,1) , uxyf(i,1)
            endif
        enddo
    enddo
    close(177)	
1312 format(9(1x,e16.8))
	 
    return
    end