program lsq

!======================================================================    
! script is used to calculate the Amplitude using Least Square Method
! data: "ondes.txt", "opa2015a.eops"
!
! work in progress
!
!======================================================================    

  use arg_nut
  use mod_matrix
  implicit none

!======================================================================  
! Variables declaration
!======================================================================  

  integer,parameter::Nbr_of_point=5980,Nbr_of_parameter=42
  integer, dimension(Nbr_of_parameter) :: l,ls,F,D,Om !integer multiplicativ coefficent to process modeling pulsation
  real (kind=xi), dimension(Nbr_of_parameter) :: ReREN,ImREN,ReMHB,ImMHB, &
    sigma,phi
  real(kind=xi) :: Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa, lpar1,lpar2,lpar3,lspar1, &
    lspar2, lspar3,Fpar1,Fpar2,Fpar3,Dpar1,Dpar2,Dpar3,omegapar1, omegapar2, &
    omegapar3, ReCOR,ImCOR,ReADD,ImADD,ReRET,ImRET,Periode
  character(10) :: Aj

  real(kind=xi), dimension(Nbr_of_point) :: t, dX, dY, errdX, errdY, corrdXdY
  real(kind=xi), dimension(2*Nbr_of_point) :: dXdY
  real(kind=xi) :: var
  character(20) :: carc

  real(kind=xi), dimension(2*Nbr_of_point,2*Nbr_of_parameter) :: M
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_parameter) :: MM, P
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_point) :: Q, MMM
  real(kind=xi), dimension(2*Nbr_of_parameter) :: Ampl
  real(kind=xi), dimension(Nbr_of_parameter) :: A, B
  integer :: i, j, k, r, s

!Read "ondes.txt" data  
  open (unit=15, file="data/ondes.txt", status='old')
  read(15,*)
  read(15,*)
  do i=1,Nbr_of_parameter
    read (15,*) Aj,l(i),ls(i),F(i),D(i),Om(i),Me,Ve,Te,Ma,Ju, &
    Sa,Ur,Ne,Pa,Periode,ReREN(i),ImREN(i),ReMHB(i),ImMHB(i),ReCOR, &
    ImCOR,ReADD,ImADD,ReRET,ImRET
!Calculate the total frequency and phase      
    sigma(i) = l(i)*sigmal + ls(i)*sigmals + F(i)*sigmaF + D(i)*sigmaD + &
            Om(i)*sigmaomega
    phi(i) = l(i)*phil + ls(i)*phils + F(i)*phiF + D(i)*phiD + Om(i)*phiomega
  end do
  close(15)

!Read "opa2015a.eops" data  
  open (unit=10, file="data/opa2015a.eops", status='old')
  do k=1,38
    read(10,*)
  end do
  do j=1,Nbr_of_point
    read (10,*) t(j), var, var, var, dX(j), dY(j), var, var, &
         carc, errdX(j), errdY(j), carc, var, var, var, corrdXdY(j), &
         var, carc, var, var, var, var, var, var, var, var, &
         var, var, var, var, carc
    if (errdX(j) .lt. 0.0009) then
       errdX(j) = 0.0001 !FIXME why 0.00001 ? and not 10000 ???
    end if
    if (errdY(j) .lt. 0.0009) then
       errdY(j) = 0.0001 !FIXME same thing 
    end if
    s = Nbr_of_point + j
    dXdY(j) = dX(j)
    dXdY(s) = dY(j)
  end do
  close(10)

!  print*, size(dX,1)
!  print*, 'Sigma (1) : ', sigma(1)
!  print*, 'Phi (1) : ', phi(1)
!  print*, 't (1) :', t(1)
!  print*, 'Uncertainty of celestial pole offset dX (1): ', errdX(1)
!  print*, 'Uncertainty of celestial pole offset dY (1): ', errdY(1)

!======================================================================  
!Least Square Method
!======================================================================  

!Create a Matrix of parameter sigma and phi with the uncertainty of dX and dY   

  do j=1,Nbr_of_point
    do k=1,Nbr_of_parameter
      M(j,k) = (1./(errdX(j)**2))*cos(sigma(k)*t(j)+phi(k))
    end do
    do k=1,Nbr_of_parameter
      r = Nbr_of_parameter + k
      M(j,r) = -(1./(errdX(j)**2))*sin(sigma(k)*t(j)+phi(k))
    end do
  end do
  do j = 1,Nbr_of_point
    s = Nbr_of_point + j
    do k=1,Nbr_of_parameter
      M(s,k) = (1./(errdY(j)**2))*sin(sigma(k)*t(j)+phi(k))
    end do
    do k=1,Nbr_of_parameter
      r = Nbr_of_parameter + k
      M(s,r) = (1./(errdY(j)**2))*cos(sigma(k)*t(j)+phi(k))
    end do
  end do

!  print*, 'Sigma (Nbr_of_parameter) : ', sigma(Nbr_of_parameter)
!  print*, 'Phi (Nbr_of_parameter) : ', phi(Nbr_of_parameter)
!  print*, 't (Nbr_of_point) :', t(Nbr_of_point)
!  print*, 'Uncertainty of celestial pole offset dX (Nbr_of_point): ', errdX(Nbr_of_point)
!  print*, 'Uncertainty of celestial pole offset dY (Nbr_of_point): ', errdY(Nbr_of_point)
!  exa = (1/(errdX(1)**2))*cos(sigma(1)*t(1)+phi(1))
!  exb = (1/(errdY(Nbr_of_point)**2))*cos(sigma(Nbr_of_parameter)*t(Nbr_of_point)+phi(Nbr_of_parameter))
!  print*, 'exa = ', exa
!  print*, 'exb = ', exb
!  print*, 'M(1,1) = ', M(1,1)
!  print*, 'M(2*Nbr_of_point,2*Nbr_of_parameter) = ', M(2*Nbr_of_point,2*Nbr_of_parameter)

! Inverse matrix with A = ((M[t]*M)^-1)*(M[t]*X)
  Q = transpose(M) ! Transopse the Matrix
  MM = matmul(Q,M) ! Multiply the Matrix
  P = inv(MM)      ! Inverse the Matrix using function which defined
  MMM = matmul(P,Q)

! Calculate the amplitude
  Ampl = matmul(MMM,dXdY)
  do i=1,Nbr_of_parameter
    s = Nbr_of_parameter + i
    A(i) = Ampl(i)
    B(i) = Ampl(s)
  end do

end program lsq
