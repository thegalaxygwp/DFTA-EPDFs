subroutine readpur()
    use allmodules
    implicit doubleprecision (a-h,o-z)
    integer :: i,j,k,n_read,nn
    character*128 file_read
    doubleprecision,allocatable :: temp1(:)
    allocate(temp1(ni*nj))

 !   file_read='D:\Data\201912+\EMHD_linear\'
 !   open(101,file='sps_16in40x2049.dat',defaultfile=file_read)
 !   read(101,*) temp
 !   write(*,*)'reading eigen_perturbations of psi...'
 !   nn=1
 !   do j=1,16
 !       do i=1,ni
 !           ps_pur(i,j)=temp(nn)
	!		!write(*,*)ps_pur(i,j)
	!		!pause
 !           nn=nn+1
 !       enddo
	!enddo
	!
	!close(101)
	!
 !   file_read='D:\Data\201912+\EMHD_linear\'
 !   open(101,file='sbz_16in40x2049.dat',defaultfile=file_read)
 !   read(101,*) temp
 !   write(*,*)'reading eigen_perturbations of bz1...'
 !   nn=1
 !   do j=1,16
 !       do i=1,ni
 !           bz_pur(i,j)=temp(nn)
	!		!write(*,*)ps_pur(i,j)
	!		!pause
 !           nn=nn+1
 !       enddo
	!enddo
	!deallocate(temp)
	!
	!close(101)
	
	
    open(101,file='egms_ps.dat',defaultfile=trim(file_dic))
    read(101,*) temp1
    write(*,*)'reading eigen_perturbations of ps1...'
    nn=1
    do i=1,ni
        do j=1,nj
            psi_z1(i,j)=temp1(nn)
			!write(*,*)ps_pur(i,j)
			!pause
            nn=nn+1
        enddo
	enddo
    write(*,*)'read complecated'
	!pause	
	close(101)
	
    open(101,file='egms_bz.dat',defaultfile=trim(file_dic))
    read(101,*) temp1
    write(*,*)'reading eigen_perturbations of bz1...'
    nn=1
    do i=1,ni
        do j=1,nj
            mgt_z1(i,j)=temp1(nn)
			!write(*,*)ps_pur(i,j)
			!pause
            nn=nn+1
        enddo
	enddo
    write(*,*)'read complecated'
	!pause	
	close(101)
	
	deallocate(temp1)

    return
    end