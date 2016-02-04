program lsq

!======================================================================    
! script is used to calculate the Amplitude using Least Square Method
! data: "ondes.txt", "opa2015a.eops"
!
! work in progress
!
!======================================================================    

  use arg_nut

  implicit none

!======================================================================  
! Variables declaration
!======================================================================  

  integer, dimension(42) :: l,ls,F,D,Om
  real (kind=xi), dimension(42) :: ReREN,ImREN,ReMHB,ImMHB, &
    sigma,phi
  real(kind=xi) :: Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa, lpar1,lpar2,lpar3,lspar1,lspar2, &
    lspar3,Fpar1,Fpar2,Fpar3,Dpar1,Dpar2,Dpar3,omegapar1,omegapar2,omegapar3, &
    ReCOR,ImCOR,ReADD,ImADD,ReRET,ImRET,Periode
  character(10) :: Aj

  real(kind=xi), dimension(5980) :: t, dX, dY, errdX, errdY, corrdXdY
  real(kind=xi) :: var, exa, exb
  character(20) :: carc

  real(kind=xi), dimension(11960,84) :: M, P

  integer :: i, j, k, r, s

  real(kind=xi), dimension(11960) :: WORK
  real(kind=xi), dimension(11960) :: IPIV
  integer :: info


!Read "ondes.txt" data  
  open (unit=15, file="data/ondes.txt", status='old')
  read(15,*)
  read(15,*)
  do i=1,42
    read (15,*) Aj,l(i),ls(i),F(i),D(i),Om(i),Me,Ve,Te,Ma,Ju, &
    Sa,Ur,Ne,Pa,Periode,ReREN(i),ImREN(i),ReMHB(i),ImMHB(i),ReCOR, &
    ImCOR,ReADD,ImADD,ReRET,ImRET
!Calculate the total frequency and phase      
    sigma(i) = l(i)*sigmal + ls(i)*sigmals + F(i)*sigmaF + D(i)*sigmaD + Om(i)*sigmaomega
    phi(i) = l(i)*phil + ls(i)*phils + F(i)*phiF + D(i)*phiD + Om(i)*phiomega
  end do
  close(15)

!Read "opa2015a.eops" data  
  open (unit=10, file="data/opa2015a.eops", status='old')
  do k=1,38
    read(10,*)
  end do
  do j=1,5980
    read (10,*) t(j), var, var, var, var, var, var, var, &
         carc, errdX(j), errdY(j), carc, var, var, var, corrdXdY(j), &
         var, carc, var, var, var, var, var, var, var, var, &
         var, var, var, var, carc
  end do
  close(10)

!  print*, 'Sigma (1) : ', sigma(1)
!  print*, 'Phi (1) : ', phi(1)
!  print*, 't (1) :', t(1)
!  print*, 'Uncertainty of celestial pole offset dX (1): ', errdX(1)
!  print*, 'Uncertainty of celestial pole offset dY (1): ', errdY(1)

!======================================================================  
!Least Square Method
!======================================================================  

!Create a Matrix of parameter sigma and phi with the uncertainty of dX and dY   
  do j=1,5980
    do k=1,42
      M(j,k) = (1/(errdX(j)**2))*cos(sigma(k)*t(j)+phi(k))
    end do
    do k=1,42
      r = 42 + k
      M(j,r) = -(1/(errdX(j)**2))*sin(sigma(k)*t(j)+phi(k))
    end do
  end do

  do j = 1,5981
    s = 5980 + j
    do k=1,42
      M(s,k) = (1/(errdY(j)**2))*sin(sigma(k)*t(j)+phi(k))
    end do
    do k=1,42
      r = 42 + k
      M(s,r) = (1/(errdY(j)**2))*cos(sigma(k)*t(j)+phi(k))
    end do
  end do

!  print*, 'Sigma (42) : ', sigma(42)
!  print*, 'Phi (42) : ', phi(42)
!  print*, 't (5980) :', t(5980)
!  print*, 'Uncertainty of celestial pole offset dX (5980): ', errdX(5980)
!  print*, 'Uncertainty of celestial pole offset dY (5980): ', errdY(5980)

  exa = (1/(errdX(1)**2))*cos(sigma(1)*t(1)+phi(1))
  exb = (1/(errdY(5980)**2))*cos(sigma(42)*t(5980)+phi(42))
  
!  print*, 'exa = ', exa
!  print*, 'exb = ', exb

!  print*, 'M(1,1) = ', M(1,1)
!  print*, 'M(11960,84) = ', M(11960,84)

!Create an inverse Matrix   
!  P = inv(M)
  
  contains 

  function inv(A) result(Ainv) ! Define the function for inversing the Matrix

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(size(A,1),size(A,2)) :: Ainv

    double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

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

end program lsq
