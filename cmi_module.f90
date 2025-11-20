	module cmi_module
	use allmodules
	use omp_lib
	use mkl_service


	contains
	!=============================================================================================!
	subroutine complex_matrix_inverse_omp(A, n_A)
	implicit none
	doublecomplex, allocatable :: A(:,:), invA(:,:), LU(:,:), L(:,:), U(:,:)
	doublecomplex, allocatable :: work(:), e(:)
	doubleprecision :: start_time, end_time
	integer :: n_A, i, j, k, error_flag
	doublecomplex :: sum_val
	allocate(invA(n_A,n_A), LU(n_A,n_A), L(n_A,n_A), U(n_A,n_A), work(n_A), e(n_A))
	LU = A  ! 复制原始矩阵

	call lu_decomposition(LU, n_A, work, error_flag)
	if (error_flag /= 0) then
		print *, "Error: Matrix is singular!"
		stop
	end if

	! 步骤2: 使用LU分解并行求解逆矩阵
	!$OMP PARALLEL DO PRIVATE(e, j, k, sum_val) SHARED(invA, LU) NUM_THREADS(max_thread)
	do i = 1, n_A
		! 初始化标准单位向量
		e = cmplx(0.0d0, 0.0d0)
		e(i) = cmplx(1.0d0, 0.0d0)

		! 前向替换 (解 Ly = e)
		do j = 1, n_A
			sum_val = e(j)
			do k = 1, j - 1
				sum_val = sum_val - LU(j, k) * e(k)
			end do
			e(j) = sum_val
		end do

		! 后向替换 (解 Ux = y)
		do j = n_A, 1, -1
			sum_val = e(j)
			do k = j + 1, n_A
				sum_val = sum_val - LU(j, k) * invA(k, i)
			end do
			invA(j, i) = sum_val / LU(j, j)
		end do
	end do
	!$OMP END PARALLEL DO
	A = invA

	deallocate(invA, LU, L, U, work, e)
	end subroutine complex_matrix_inverse_omp

	subroutine complex_matrix_inverse_omp_di(LU, n_A)
	implicit none
	doublecomplex, allocatable :: invA(:,:), LU(:,:)
	doublecomplex, allocatable :: e(:)
	doubleprecision :: start_time, end_time
	integer :: n_A, i, j, k
	doublecomplex :: sum_val
	allocate(invA(n_A,n_A), e(n_A))
	e = cmplx(0.0d0, 0.0d0)

	! 步骤1: 使用LU分解并行求解逆矩阵
	!$omp parallel shared(n_A) num_threads(max_thread)
	!$omp do collapse(1) private(e, j, k, sum_val)! schedule(dynamic,Nr2_dy)
	do i = 1, n_A
		! 初始化标准单位向量
		e(i) = cmplx(1.0d0, 0.0d0)
		! 前向替换 (解 Ly = e)
		do j = 1, n_A
			sum_val = e(j)
			do k = 1, j - 1
				sum_val = sum_val - LU(j, k) * e(k)
			end do
			e(j) = sum_val
		end do
		! 后向替换 (解 Ux = y)
		do j = n_A, 1, -1
			sum_val = e(j)
			do k = j + 1, n_A
				sum_val = sum_val - LU(j, k) * invA(k, i)
			end do
			invA(j, i) = sum_val / LU(j, j)
		end do
	end do
	!$omp enddo
	!$omp end parallel


	!!$OMP PARALLEL DO PRIVATE(e, j, k, sum_val) SHARED(invA, LU) NUM_THREADS(max_thread)
	!do i = 1, n_A
	!    ! 初始化标准单位向量
	!    e = cmplx(0.0d0, 0.0d0)
	!    e(i) = cmplx(1.0d0, 0.0d0)
	!    ! 前向替换 (解 Ly = e)
	!    do j = 1, n_A
	!        sum_val = e(j)
	!        do k = 1, j - 1
	!            sum_val = sum_val - LU(j, k) * e(k)
	!        end do
	!        e(j) = sum_val
	!    end do
	!    ! 后向替换 (解 Ux = y)
	!    do j = n_A, 1, -1
	!        sum_val = e(j)
	!        do k = j + 1, n_A
	!            sum_val = sum_val - LU(j, k) * invA(k, i)
	!        end do
	!        invA(j, i) = sum_val / LU(j, j)
	!    end do
	!end do
	!!$OMP END PARALLEL DO
	LU = invA

	deallocate(invA, e)
	end subroutine complex_matrix_inverse_omp_di

	! LU分解 (带部分主元选取)
	subroutine lu_decomposition(mat, n_A, work, error_flag)
	doublecomplex, intent(inout) :: mat(:,:)
	integer, intent(in) :: n_A
	doublecomplex, intent(out) :: work(:)
	integer, intent(out) :: error_flag
	integer :: i, j, k, pivot_row
	doublecomplex :: temp, pivot
	doubleprecision :: max_val

	error_flag = 0
	do j = 1, n_A
		! 查找主元
		pivot_row = j
		max_val = abs(mat(j,j))
		do i = j + 1, n_A
			if (abs(mat(i,j)) > max_val) then
				max_val = abs(mat(i,j))
				pivot_row = i
			end if
		end do

		! 检查奇异性
		if (abs(mat(pivot_row,j)) < 1e-12) then
			error_flag = j
			return
		end if

		! 行交换
		if (pivot_row /= j) then
			work = mat(j,:)
			mat(j,:) = mat(pivot_row,:)
			mat(pivot_row,:) = work
		end if

		! 高斯消元
		pivot = mat(j,j)
		do i = j + 1, n_A
			mat(i,j) = mat(i,j) / pivot
			mat(i,j+1:n_A) = mat(i,j+1:n_A) - mat(i,j) * mat(j,j+1:n_A)
		end do
	end do
	end subroutine lu_decomposition


	!  =====================================================================

	SUBROUTINE zgetri_hand( N, A, LDA, IPIV, WORK, LWORK, INFO )
	!	  LAPACK computational routine --
	!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
	!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	integer :: INFO, LDA, LWORK, N
	integer :: IPIV( * )
	doublecomplex ::  A( LDA, * ), WORK( * )
	doublecomplex :: ZERO, ONE
	LOGICAL :: LQUERY
	integer :: I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,NBMIN, NN
	integer :: ILAENV

	zero = ( 0.0d+0, 0.0d+0 )
	one = ( 1.0d+0, 0.0d+0 )
	info = 0
	nb = ilaenv( 1, 'ZGETRI', ' ', n, -1, -1, -1 )
	lwkopt = max( 1, n*nb )
	work( 1 ) = lwkopt
	lquery = ( lwork.EQ.-1 )
	IF( n.LT.0 ) THEN
		info = -1
	ELSE IF( lda.LT.max( 1, n ) ) THEN
		info = -3
	ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
		info = -6
	END IF

	IF( info.NE.0 ) THEN
		CALL xerbla( 'ZGETRI', -info )
		RETURN
	ELSE IF( lquery ) THEN
		RETURN
	END IF

	IF( n.EQ.0 )  RETURN

	CALL ztrtri( 'Upper', 'Non-unit', n, a, lda, info )
	IF( info.GT.0 )   RETURN

	nbmin = 2
	ldwork = n
	IF( nb.GT.1 .AND. nb.LT.n ) THEN
		iws = max( ldwork*nb, 1 )
		IF( lwork.LT.iws ) THEN
			nb = lwork / ldwork
			nbmin = max( 2, ilaenv( 2, 'ZGETRI', ' ', n, -1, -1,-1 ) )
		END IF
	ELSE
		iws = n
	END IF

	IF( nb.LT.nbmin .OR. nb.GE.n ) THEN
		DO 20 j = n, 1, -1
			DO 10 i = j + 1, n
				work( i ) = a( i, j )
				a( i, j ) = zero
10			CONTINUE
			IF( j.LT.n )   CALL zgemv( 'No transpose', n, n-j, -one, a( 1, j+1 ),lda, work( j+1 ), 1, one, a( 1, j ), 1 )
20		CONTINUE
	ELSE
		nn = ( ( n-1 ) / nb )*nb + 1
		DO 50 j = nn, 1, -nb
			jb = min( nb, n-j+1 )
			DO 40 jj = j, j + jb - 1
				DO 30 i = jj + 1, n
					work( i+( jj-j )*ldwork ) = a( i, jj )
					a( i, jj ) = zero
30				CONTINUE
40			CONTINUE
			IF( j+jb.LE.n )   CALL zgemm( 'No transpose', 'No transpose', n, jb,n-j-jb+1, -one, a( 1, j+jb ), lda,work( j+jb ), ldwork, one, a( 1, j ), lda )
			CALL ztrsm( 'Right', 'Lower', 'No transpose', 'Unit', n,jb,one, work( j ), ldwork, a( 1, j ), lda )
50		CONTINUE
	END IF
	DO 60 j = n - 1, 1, -1
		jp = ipiv( j )
		IF( jp.NE.j )   CALL zswap( n, a( 1, j ), 1, a( 1, jp ), 1 )
60	CONTINUE
	work( 1 ) = iws
	RETURN

	END SUBROUTINE zgetri_hand
	
	
	!  =====================================================================

      SUBROUTINE zgetrf_hand( M, N, A, LDA, IPIV, INFO )
      integer ::INFO, LDA, M, N
      integer ::IPIV( * )
      doublecomplex :: A( LDA, * )
      doublecomplex :: ONE
      integer ::I, IINFO, J, JB, NB
      integer ::ILAENV
	  
	  one = ( 1.0d+0, 0.0d+0 )
	  
	  info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGETRF', -info )
         RETURN
      END IF
      IF( m.EQ.0 .OR. n.EQ.0 )   RETURN
      nb = ilaenv( 1, 'ZGETRF', ' ', m, n, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.min( m, n ) ) THEN
         CALL zgetrf2( m, n, a, lda, ipiv, info )
      ELSE
         DO 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )
            CALL zgetrf2( m-j+1, jb, a( j, j ), lda, ipiv( j ),iinfo )
            IF( info.EQ.0 .AND. iinfo.GT.0 )        info = iinfo + j - 1
            DO 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       CONTINUE
            CALL zlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
            IF( j+jb.LE.n ) THEN
               CALL zlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,ipiv, 1 )
               CALL ztrsm( 'Left', 'Lower', 'No transpose', 'Unit',jb,n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),lda )
               IF( j+jb.LE.m ) THEN
                  CALL zgemm( 'No transpose', 'No transpose',m-j-jb+1,n-j-jb+1, jb, -one, a( j+jb, j ), lda,a( j, j+jb ), lda, one, a( j+jb, j+jb ),lda )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN

	  END SUBROUTINE zgetrf_hand
	  
	  
	end module cmi_module