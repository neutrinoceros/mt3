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
  use mod_series
  implicit none

!======================================================================  
! Variables declaration
!======================================================================  

  integer,parameter::Nbr_of_point=5980,Nbr_of_parameter=42
  integer, dimension(Nbr_of_parameter) :: l,ls,F,D,Om                          !integer multiplicative coefficent to process modeling pulsation
  real (kind=xi), dimension(Nbr_of_parameter) :: ReREN,ImREN,ReMHB,ImMHB,    &
    sigma,phi                                                                  !reading variables not used in the code
  real(kind=xi) :: Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa, lpar1,lpar2,lpar3,lspar1,     &
    lspar2, lspar3,Fpar1,Fpar2,Fpar3,Dpar1,Dpar2,Dpar3,omegapar1, omegapar2, &
    omegapar3, ReCOR,ImCOR,ReADD,ImADD,ReRET,ImRET,Periode                     !reading variables not used in the code
  character(10) :: Aj                                                          !reading variables not used in the code

  real(kind=xi), dimension(Nbr_of_point)   :: t              !time of the mesure in julian days
  real(kind=xi), dimension(Nbr_of_point)   :: dX, dY         !delta between observe nutation and IAU2000 in mas
  real(kind=xi), dimension(Nbr_of_point)   :: errdX, errdY   !error of the mesurment in mas
  real(kind=xi), dimension(Nbr_of_point)   :: corrdXdY       !correlation of the measure
  real(kind=xi), dimension(2*Nbr_of_point) :: dXdY           !matrix of observable
  real(kind=xi) :: var                                       !reading variable
  character(20) :: carc                                      !reading variable



  real(kind=xi), dimension(2*Nbr_of_point,2*Nbr_of_parameter + 4) :: M             !matrix like dXdY = M Ampl
  real(kind=xi), dimension(2*Nbr_of_parameter + 4,2*Nbr_of_parameter + 4) :: MM, P !MM is equal to tranpose(M) matrix product M
  real(kind=xi), dimension(2*Nbr_of_parameter + 4,2*Nbr_of_point) :: Q, MMM
  real(kind=xi), dimension(2*Nbr_of_parameter + 4) :: Ampl                         !complex amplitude that we need to adjust
  real(kind=xi), dimension(2*Nbr_of_point) :: M_dot_amp, tmpMa !matrix to compute the error
  real(kind=xi) ::sigma_lsq ! residu


  real(kind=xi), dimension(Nbr_of_parameter) :: A, B
  real(kind=xi), dimension(Nbr_of_parameter) :: errA,errB
  integer :: i, j, k, r, s                                                 !iterators
  real(kind=xi) :: phase                                                   !instant phase used in matrix computation loop
  complex(kind=xi), dimension(Nbr_of_point)   :: SrA, SrB, Sr              !results of the forier series of the amplitudes (series and MHB)

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
!This matrix is composed of 4 distinct cadrans

  do j=1,Nbr_of_point
     s = Nbr_of_point + j
     do k=1,Nbr_of_parameter
        r = Nbr_of_parameter + k
        phase  = sigma(k)*t(j)+phi(k)
        M(j,k) = +1./errdX(j)*cos(phase)
        M(j,r) = -1./errdX(j)*sin(phase)
        M(s,k) = +1./errdY(j)*sin(phase)
        M(s,r) = +1./errdY(j)*cos(phase)
     end do
  end do


  ! addition of more parameters : systematic slope and bias

  M(1              : Nbr_of_point   , 2*Nbr_of_parameter + 1) = t
  M(1              : Nbr_of_point   , 2*Nbr_of_parameter + 2) = 1.0_xi

  M(Nbr_of_point+1 : 2*Nbr_of_point , 2*Nbr_of_parameter + 3) = t
  M(Nbr_of_point+1 : 2*Nbr_of_point , 2*Nbr_of_parameter + 4) = 1.0_xi



! Inverse matrix with A = ((M[t]*M)^-1)*(M[t]*X)
  Q = transpose(M) ! transpostion
  MM = matmul(Q,M) ! matrix product
  P = inv(MM)      ! matrix inversion using function defined in mod_matrix.f90
  MMM = matmul(P,Q)


! Calculate corrections to complex amplitudes in the model
  Ampl = matmul(MMM,dXdY)

  ! processing the error on lsq
  M_dot_amp = matmul(M,Ampl)
  tmpMa = (dXdY-M_dot_amp)**2
  sigma_lsq = sqrt(sum(tmpMa))/(Nbr_of_point-Nbr_of_parameter)

  open (unit=12,file="amplitude.dat",status="replace")
  write(12,fmt='(5 A26)') "#dX", "dY", "dX_error", "dY_error", "Frequency"
  do i=1,Nbr_of_parameter
    s = Nbr_of_parameter + i

    A(i)    = Ampl(i)
    B(i)    = Ampl(s)
    errA(i) = sqrt(P(i,i)) * sigma_lsq
    errB(i) = sqrt(P(s,s)) * sigma_lsq

    write(12,fmt='(5 E26.16)') A(i), B(i), errA(i), errB(i),sigma(i)
  end do
  close(unit=12)

! Calculate the value of the X and Y from series and MHB 
  open (unit=13,file="series.dat",status="replace")
  write(13,fmt='(5 A26)') "#dX (series + MHB)", "dY (series + MHB)", "dX (observation)", "dY (observation)", "time"
  do j = 1, Nbr_of_point
    SrA(j) = ser(Nbr_of_parameter,A,B,sigma,phi,t(j))         ! X and Y for the series
    SrB(j) = ser(Nbr_of_parameter,ReREN,ImREN,sigma,phi,t(j)) ! X and Y for the MHB
    Sr(j)  = SrA(j) !+ SrB(j)                                  ! the sum of X and Y from series and MHB
    write(13,fmt='(5 E26.16)') real(Sr(j)), aimag(Sr(j)), dX(j), dY(j), t(j)
  end do
  close(unit=13)

  ! printing of slope and bias fitting results
  s = 2*Nbr_of_parameter + 1
  print*,
  print*, 'slope (a), and bias (b) :  (a_re, b_re, a_im, b_im) :'
  print*, Ampl(s:s+3)

  print*,
end program lsq
