	program DFTA
    use allmodules
    use dfport
	use function_math
    character*128 infile1,infile2,infile3
    logical istatus1,istatus2
    
    tstart=omp_get_wtime()
    !$omp parallel !single
    num_core    =OMP_GET_NUM_PROCS()
    max_thread  =OMP_GET_MAX_THREADS()
    !$omp end parallel !single
    !pause
    WRITE(*,*)'we have',num_core,'cores, and',max_thread,'threads can be usesd!'

    !call omp_set_nested(1)
	call omp_set_max_active_levels(2)
	call omp_set_dynamic(1)
	
    write(infile2,'(a27)')'E:\585\2025\emhddata\202510'  !   Set Output File
    filedata=trim(infile2)
    write(*,*)filedata
    infile1='MD '//trim(infile2)    !   The statement to create a folder uses the md command, first defining it as a string.
    infile3='CD '//trim(infile2)    !   To check if a folder exists using the cd command:
    istatus1=system(infile3)        !   Check if folder exists
    if(istatus1)istatus2=system(infile1)
	
	call transprt()        
    call dimens_allocate()          !   Allocate The Dimens
    call set_xy()                   !   Set Operator DFT-D Matrix Of 1-3 Order
    call poutparam()				!   Output The Caculation Parameters	!
    call equation_set()             !	Set Initial Function	!
    call big_matrix_set()			!	Set A Function	!
	call analysies_DFTA()           !   DFTA Caculate Test And Report Results
		
    tend=omp_get_wtime()
    write(*,*)tend-tstart
    pause
    end program DFTA