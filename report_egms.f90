    subroutine report_egms()
    use allmodules
    implicit doubleprecision (a-h,o-z)
    integer i,j
    character*15 fdata
    !=============================================================================================!
    open(311,file='egms_ps.dat',defaultfile=trim(file_dic))
    open(312,file='egms_bz.dat',defaultfile=trim(file_dic))
    !=============================================================================================!
    do i=1,ni
        do j=1,nj
			write(311,1312)1e-6*psi_z1(i,j)
			write(312,1312)1e-6*mgt_z1(i,j)
        enddo
    enddo
    !=============================================================================================!
    close(311)
    close(312)
1312 format(1(e18.6e6))
    return
    end