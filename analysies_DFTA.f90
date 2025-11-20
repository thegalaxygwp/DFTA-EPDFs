    subroutine analysies_DFTA()
    use module_solvers
    implicit doubleprecision (a-h,o-z)
    integer :: i,j,k	
    character*128 outfile1,outfile2
    !=============================================================================================!
    if(n_ka/=0)then
        if(n_reduce==0)then
            !call trans_sin(fxyf_DFT,fxyf)
            !call trans_cos(hxyf_DFT,hxyf)
            call solve_x_fast_linear(gxyf_DFT,fxyf)
            call solve_y_fast_linear(uxyf_DFT,hxyf)
        elseif(n_reduce/=0)then
            call solve_x(gxyf_DFT,fxyf)
            call solve_y(uxyf_DFT,hxyf)
            !call trans_sin(fxyf_DFT,fxyf)
            !call trans_cos(hxyf_DFT,hxyf)
        endif
    elseif(kappa==0 .and. n_ka==0)then
        call trans_sin(fxyf_DFT,fxyf)
        call trans_cos(hxyf_DFT,hxyf)
        call solve_uniform_sin(gxyf_DFT,fxyf_DFT)
        call solve_uniform_cos(uxyf_DFT,hxyf_DFT)
	endif
	
	call inves_sin(gxyf_DFT,gxyf)
    call inves_cos(uxyf_DFT,uxyf)
	
	write(outfile1,'(a13,i4,a1,i2,a1,i1,a4)')'function_data-',ni,'x',nj,'_',n_reduce,'.dat'  !   Set Output File
    !write(outfile2,'(a13,i4,a1,i2,a1,i1,a4)')'function_dfta-',ni,'x',nj,'_',n_reduce,'.dat'  !   Set Output File

	open(177,file=outfile1,defaultfile=trim(filedata))
    write(177,*)'x y gxyf0 uxyf0 fxyf hxyf sxyf gxyf uxyf dgxy duxy'
    do i=1,ni
        do j=1,nj+1
            if(j<=nj)then
                write(177,1312) cod_x(i,j),cod_y(i,j), &
                    gxyf0(i,j),uxyf0(i,j), &
                    fxyf(i,j) ,hxyf(i,j) , &
					sxyf(i,j) ,gxyf(i,j) , uxyf(i,j), gxyf0(i,j)-gxyf(i,j), uxyf0(i,j)-uxyf(i,j)
            elseif(j==nj+1)then
                write(177,1312) cod_x(i,1),cod_y(i,1), &
                    gxyf0(i,1),uxyf0(i,1), &
                    fxyf(i,1) ,hxyf(i,1) , &
					sxyf(i,1) ,gxyf(i,1) , uxyf(i,1), gxyf0(i,1)-gxyf(i,1), uxyf0(i,1)-uxyf(i,1)
            endif
        enddo
    enddo
    close(177)	
1312 format(11(1x,e16.8))


!	open(177,file=outfile2,defaultfile=trim(filedata))
!    !write(177,*)'x y gxyf_DFT gxyf_DFT uxyf_DFT uxyf_DFT'
!    do i=1,nr
!        do j=1,ns
!            write(177,1313) i,j,real(gxyf_DFT(i,j)),aimag(gxyf_DFT(i,j)),real(uxyf_DFT(i,j)),aimag(uxyf_DFT(i,j))
!        enddo
!    enddo
!    close(177)	
!1313 format(6(1x,e16.8))
    !=============================================================================================!
    return
    end