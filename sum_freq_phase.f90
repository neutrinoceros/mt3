module sum_freq_phase

!======================================================================    
! script used to calculate the sum of frequency and phase 
! using the data in "ondes.txt"
!======================================================================    

  implicit none

!======================================================================  
! Variables declaration
!======================================================================  

  integer, parameter :: xi = selected_real_kind (15)
  integer, dimension(:), allocatable :: l,ls,F,D,Om,Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa 
  real (kind=xi), dimension(:), allocatable :: Periode,ReREN, ImREN,ReMHB,ImMHB,ReCOR, &
    ImCOR,ReADD,ImADD,ReRET,ImRET, sigma, phi
  real(kind=xi) :: phil, phils, phiF, phiD, phiomega, &
    sigmal, sigmals, sigmaF, sigmaD, sigmaomega, lpar1, lpar2, lpar3, lspar1, lspar2, &
    lspar3, Fpar1, Fpar2, Fpar3, Dpar1, Dpar2, Dpar3, omegapar1, omegapar2, omegapar3
  character(10) :: Aj
  integer :: i, line=42
  
  real (kind=xi), parameter :: pi=4._xi*atan(1._xi)
   
  allocate(l(line))
  allocate(ls(line))
  allocate(F(line))
  allocate(D(line))
  allocate(Om(line))
  allocate(Me(line))
  allocate(Ve(line))
  allocate(Te(line))
  allocate(Ma(line))
  allocate(Ju(line))
  allocate(Sa(line))
  allocate(Ur(line))
  allocate(Ne(line))
  allocate(Pa(line))
  allocate(Periode(line))
  allocate(ReREN(line))
  allocate(ImREN(line))
  allocate(ReMHB(line))
  allocate(ImMHB(line))
  allocate(ReCOR(line))
  allocate(ImCOR(line))
  allocate(ReADD(line))
  allocate(ImADD(line))
  allocate(ReRET(line))
  allocate(ImRET(line))
  allocate(sigma(line))
  allocate(phi(line))

  !Define the parameter of Fundamental Argument of Nutation (in radian)
  !phi.. (phase) equal to constant part of the equation
  phil = 134.96340251_xi*pi/180._xi
  phils = (357.52910918_xi*pi/180._xi)
  phiF = 93.27209062_xi*pi/180._xi
  phiD = (297.85019547_xi*pi/180._xi)
  phiomega = 125.04455501_xi*pi/180._xi

  !sigma.. (frequency) equal to coeficient of "t" of the equation
  sigmal =(1717915923.2178_xi/3600._xi)*pi/180._xi
  sigmals =(129596581.0481_xi/3600._xi)*pi/180._xi
  sigmaF =(1739527262.8478_xi/3600._xi)*pi/180._xi
  sigmaD =(1602961601.2090_xi/3600._xi)*pi/180._xi
  sigmaomega =(-6962890.5431_xi/3600._xi)*pi/180._xi

  !Read "ondes.txt" data  
  open (unit=15, file="ondes.txt", status='old')
  read(15,*)
  read(15,*)
  do i=1,42
    read (15,*) Aj,l(i),ls(i),F(i),D(i),Om(i),Me(i),Ve(i),Te(i),Ma(i),Ju(i), &
    Sa(i),Ur(i),Ne(i),Pa(i),Periode(i),ReREN(i),ImREN(i),ReMHB(i),ImMHB(i),ReCOR(i), &
    ImCOR(i),ReADD(i),ImADD(i),ReRET(i),ImRET(i)

  !Calculate the total frequency and phase      
    sigma(i) = l(i)*sigmal + ls(i)*sigmals + F(i)*sigmaF + D(i)*sigmaD + Om(i)*sigmaomega
    phi(i) = l(i)*phil + ls(i)*phils + F(i)*phiF + D(i)*phiD + Om(i)*phiomega

  end do
  close(15)

  deallocate(l)
  deallocate(ls)
  deallocate(F)
  deallocate(D)
  deallocate(Om)
  deallocate(Me)
  deallocate(Ve)
  deallocate(Te)
  deallocate(Ma)
  deallocate(Ju)
  deallocate(Sa)
  deallocate(Ur)
  deallocate(Ne)
  deallocate(Pa)
  deallocate(Periode)
  deallocate(ReREN)
  deallocate(ImREN)
  deallocate(ReMHB)
  deallocate(ImMHB)
  deallocate(ReCOR)
  deallocate(ImCOR)
  deallocate(ReADD)
  deallocate(ImADD)
  deallocate(ReRET)
  deallocate(ImRET)
  deallocate(sigma)
  deallocate(phi)

end module sum_freq_phase
