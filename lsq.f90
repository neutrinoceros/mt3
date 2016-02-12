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
    sigma,phi !reading variable not use in the code
  real(kind=xi) :: Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa, lpar1,lpar2,lpar3,lspar1, &
    lspar2, lspar3,Fpar1,Fpar2,Fpar3,Dpar1,Dpar2,Dpar3,omegapar1, omegapar2, &
    omegapar3, ReCOR,ImCOR,ReADD,ImADD,ReRET,ImRET,Periode !reading variable not use in the code
  character(10) :: Aj !reading variable not use in the code

  real(kind=xi), dimension(Nbr_of_point) :: t !time of the mesure in julian days
  real(kind=xi), dimension(Nbr_of_point) :: dX, dY !delta between observe nutation and IAU2000 in mas
  real(kind=xi), dimension(Nbr_of_point) :: errdX, errdY !error of the mesurment in mas
  real(kind=xi), dimension(Nbr_of_point) :: corrdXdY !correlation of the measure
  real(kind=xi), dimension(2*Nbr_of_point) :: dXdY !matrix of observable
  real(kind=xi) :: var !reading variable
  character(20) :: carc !reading variable

  real(kind=xi), dimension(2*Nbr_of_point,2*Nbr_of_parameter) :: M !matrix like dXdY = M Ampl
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_parameter) :: MM, P !MM is equal to tranpose(M) matrix product M
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_point) :: Q, MMM
  real(kind=xi), dimension(2*Nbr_of_parameter) :: Ampl !complex amplitude that we need to adjust
  real(kind=xi), dimension(Nbr_of_parameter) :: A, B
  integer :: i, j, k, r, s !loop variable

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
    
    !some data have 0.000 as error that meen tha the value is not true
    !we change this error in 100000. to kill this point
    if (errdX(j) .lt. 0.0009) then
       errdX(j) = 1.e6_xi
    end if
    if (errdY(j) .lt. 0.0009) then
       errdY(j) = 1.e6_xi
    end if

  end do
  close(10)

  !computation of the observationnal array
  dXdY(:Nbr_of_point)=dX/errdX
  dXdY(Nbr_of_point+1:)=dY/errdY

  !conversion of time array : days -->julian century
  t=(t-t2000)/Nbr_days_in_Century

  ! do i=1,Nbr_of_point,1
  !   print*,t(i)*36525._xi
  ! end do

!======================================================================  
!Least Square Method
!======================================================================  

!Create a Matrix of parameter sigma and phi with the uncertainty of dX and dY   
!This matrix is compose by 4 distinct cadran

  do j=1,Nbr_of_point
    !processing the 1st cadran of the matrix
    do k=1,Nbr_of_parameter
      M(j,k) = 1./errdX(j)*cos(sigma(k)*t(j)+phi(k))
    end do

    !processing the 2nd cadran of the matrix
    do k=1,Nbr_of_parameter
      r = Nbr_of_parameter + k
      M(j,r) = -1./errdX(j)*sin(sigma(k)*t(j)+phi(k))
    end do

  end do

  do j = 1,Nbr_of_point
    s = Nbr_of_point + j

    !processing the 3rd cadran of the matrix
    do k=1,Nbr_of_parameter
      M(s,k) = 1./errdY(j)*sin(sigma(k)*t(j)+phi(k))
    end do

    !processing the 4th cadran of the matrix
    do k=1,Nbr_of_parameter
      r = Nbr_of_parameter + k
      M(s,r) = 1./errdY(j)*cos(sigma(k)*t(j)+phi(k))
    end do

  end do

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
    ! print*,sigma(i)
    print*,i, A(i), B(i)
  end do

end program lsq
