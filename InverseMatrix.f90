program InverseMatrix

  implicit none

  double precision, dimension (4,4) :: C, B
  double precision, dimension (2,4) :: D
  double precision, dimension (4,4) :: E, F
  integer :: i, j

  C = 0.0

  do i=1,4
    C(i,i) = 1.0/i
  end do
  
  print*, 'Matrix :'
  do i=1,4
    print*, C(i,:)
  end do


  F = inv(C)
  print*, 'Inverse Matrix :'
  do i=1,4
    print*, F(i,:)
  end do

  contains 

! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
  function inv(A) result(Ainv)

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(size(A,1),size(A,2)) :: Ainv

    double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info


    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

end program InverseMatrix
