program complex_matrix_inverse_omp
  use omp_lib
  implicit none
  integer, parameter :: dp = kind(0.d0)  ! 双精度
  complex(dp), allocatable :: A(:,:), invA(:,:), LU(:,:), L(:,:), U(:,:)
  complex(dp), allocatable :: work(:), e(:)
  real(dp) :: start_time, end_time
  integer :: n_A, i, j, k, num_threads, error_flag
  complex(dp) :: sum_val

  ! 矩阵维度设置
  n_A = 4  ! 可修改为任意维度
  num_threads = 4  ! 设置OpenMP线程数
  
  ! 创建并初始化一个随机Hermitian复数矩阵 (保证可逆)
  allocate(A(n_A,n_A), invA(n_A,n_A), LU(n_A,n_A), L(n_A,n_A), U(n_A,n_A), work(n_A), e(n_A))
  call random_seed()
  call random_matrix(A, n_A)
  
  ! 创建Hermitian矩阵: A = A + A^H
  A = 0.5_dp * (A + conjg(transpose(A)))
  
  ! 增加对角优势确保可逆
  do i = 1, n_A
    A(i,i) = A(i,i) + cmplx(n_A, n_A, dp)
  end do

  ! 打印原始矩阵
  print *, "Original Matrix:"
  do i = 1, n_A
    print "(\*(a, f7.4, a, f7.4, a))", ('(', real(A(i,j)), ',', aimag(A(i,j)), ')  ', j = 1, n_A)
  end do

  ! 记录开始时间
  start_time = omp_get_wtime()
  
  ! 步骤1: 执行LU分解 (带部分主元选取)
  LU = A  ! 复制原始矩阵
  call lu_decomposition(LU, n_A, work, error_flag)
  
  if (error_flag /= 0) then
    print *, "Error: Matrix is singular!"
    stop
  end if

  ! 步骤2: 使用LU分解并行求解逆矩阵
  !$OMP PARALLEL DO PRIVATE(e, j, k, sum_val) SHARED(invA, LU) NUM_THREADS(num_threads)
  do i = 1, n_A
    ! 初始化标准单位向量
    e = cmplx(0.0_dp, 0.0_dp, dp)
    e(i) = cmplx(1.0_dp, 0.0_dp, dp)
    
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

  ! 记录结束时间
  end_time = omp_get_wtime()
  
  ! 打印逆矩阵
  print *, "Inverse Matrix:"
  do i = 1, n_A
    print "(\*(a, f7.4, a, f7.4, a))", ('(', real(invA(i,j)), ',', aimag(invA(i,j)), ')  ', j = 1, n_A)
  end do
  
  ! 验证结果: 计算 A * A^{-1} 应接近单位矩阵
  print *, "Verification (A \* A^{-1}):"
  call verify_inverse(A, invA, n_A)
  
  ! 打印计算时间
  print '(a, f9.6, a)', "Computation time: ", end_time - start_time, " seconds"
  print '(a, i2)', "Number of threads used: ", num_threads

  deallocate(A, invA, LU, L, U, work, e)

contains

  ! 生成随机复数矩阵
  subroutine random_matrix(mat, n_A)
    complex(dp), intent(out) :: mat(:,:)
    integer, intent(in) :: n_A
    real(dp) :: re, im
    integer :: i, j
    
    do i = 1, n_A
      do j = 1, n_A
        call random_number(re)
        call random_number(im)
        mat(i,j) = cmplx(re, im, dp)
      end do
    end do
  end subroutine random_matrix

  ! LU分解 (带部分主元选取)
  subroutine lu_decomposition(mat, n_A, work, error_flag)
    complex(dp), intent(inout) :: mat(:,:)
    integer, intent(in) :: n_A
    complex(dp), intent(out) :: work(:)
    integer, intent(out) :: error_flag
    integer :: i, j, k, pivot_row
    complex(dp) :: temp, pivot
    real(dp) :: max_val
    
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

  ! 验证逆矩阵 (计算 A * A^{-1})
  subroutine verify_inverse(A, invA, n_A)
    complex(dp), intent(in) :: A(:,:), invA(:,:)
    integer, intent(in) :: n_A
    complex(dp) :: product(n_A,n_A), diff
    real(dp) :: error
    integer :: i, j, k
    
    product = matmul(A, invA)
    
    do i = 1, n_A
      print "(\*(a, f7.4, a, f7.4, a))", ('(', real(product(i,j)), ',', aimag(product(i,j)), ')  ', j = 1, n_A)
    end do
    
    ! 计算与单位矩阵的误差
    error = 0.0_dp
    do i = 1, n_A
      do j = 1, n_A
        if (i == j) then
          diff = product(i,j) - cmplx(1.0_dp, 0.0_dp, dp)
        else
          diff = product(i,j)
        end if
        error = max(error, abs(diff))
      end do
    end do
    print '(a, e10.3)', "Max error: ", error
  end subroutine verify_inverse

end program complex_matrix_inverse_omp
