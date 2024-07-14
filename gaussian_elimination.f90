program gaussian_elimination
  implicit none
  integer :: n ! Size of the system
  real(8), allocatable :: A(:,:), b(:), x(:) ! Coefficient matrix, right-hand side vector, solution vector
  integer :: i, j, k ! `i` is the row index, `j` is the column index, `k` is the pivot index

  ! Ask for the size of the system
  print *, 'Enter the size of the system: '
  read *, n

  ! Allocate arrays
  allocate(A(n, n), b(n), x(n))

  ! Initialize matrix A and vector b in row-major format
  print *, 'Enter the coefficients of the matrix A (row by row): '
  do i = 1, n
    do j = 1, n
      read *, A(i, j)
    end do
  end do

  print *, 'Enter the right-hand side vector b: '
  do i = 1, n
    read *, b(i)
  end do

  ! Forward elimination
  do k = 1, n-1, 1
    do i = k+1, n, 1
      if (A(k, k) /= 0.0) then
        b(i) = b(i) - A(i, k) / A(k, k) * b(k)
        do j = k+1, n, 1
          A(i, j) = A(i, j) - A(i, k) / A(k, k) * A(k, j)
        end do
        A(i, k) = 0.0
      else
        print *, "Error: Division by zero"
        stop
      end if
    end do
  end do

  ! Back substitution
  x(n) = b(n) / A(n, n)
  do i = n-1, 1, -1
    x(i) = b(i)
    do j = i+1, n, 1
      x(i) = x(i) - A(i, j) * x(j)
    end do
    x(i) = x(i) / A(i, i)
  end do

  ! Output the solution
  print *, "Solution vector x: "
  do i = 1, n
    print *, x(i)
  end do

  ! Deallocate arrays
  deallocate(A, b, x)

end program gaussian_elimination
