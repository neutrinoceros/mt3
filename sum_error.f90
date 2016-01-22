program sum_error

  implicit none
  integer, parameter :: xi = selected_real_kind (15)
  integer, dimension(:), allocatable :: l,ls,F,D,Om,Me,Ve,Te,Ma,Ju,Sa,Ur,Ne,Pa 
  real (kind=xi), dimension(:), allocatable :: Periode,ReREN, ImREN,ReMHB,ImMHB,ReCOR, &
    ImCOR,ReADD,ImADD,ReRET,ImRET, sigma, phi
  real(kind=xi) :: phil, phils, phiF, phiD, phiomega, &
    sigmal, sigmals, sigmaF, sigmaD, sigmaomega, lpar1, lpar2, lspar1, lspar2, &
    Fpar1, Fpar2, Dpar1, Dpar2, omegapar1, omegapar2
  character(10) :: Aj
  integer :: i, line=42
  real (kind=xi), parameter :: pi = 3.14159265359_xi
   
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

  !Parameter of sigma
  lpar1 = 134.96340251_xi*pi/180._xi
  lspar1 = (357.52910918_xi*pi/180._xi) - 2*pi
  Fpar1 = 93.27209062_xi*pi/180._xi
  Dpar1 = (297.85019547_xi*pi/180._xi) - 2*pi
  omegapar1 = 125.04455501_xi*pi/180._xi

  lpar2 =(1717915923.2178_xi/3600._xi)*pi/180._xi
  lspar2 =(129596581.0481_xi/3600._xi)*pi/180._xi
  Fpar2 =(1739527262.8478_xi/3600._xi)*pi/180._xi
  Dpar2 =(1602961601.2090_xi/3600._xi)*pi/180._xi
  omegapar2 =(6962890.5431_xi/3600._xi)*pi/180._xi

  phil = acos(lpar1)
  print*, "1"
  sigmal = lpar2/sin(phil)
  print*, "2"
  phils = acos(lspar1)
  print*, "3"
  sigmals = lspar2/sin(phils)
  print*, "4"
  phiF = acos(Fpar1)
  print*, "5"
  sigmaF = Fpar2/sin(phiF)
  print*, "6"
  phiD = acos(Dpar1)
  print*, "7"
  sigmaD = Dpar2/sin(phiD)
  print*, "8"
  phiomega = acos(omegapar1)
  print*, "9"
  sigmaomega = omegapar2/sin(phiomega)
  print*, "10"

  print*, phil, sigmal, phils, sigmals, phiF, sigmaF, phiD, sigmaD, phiomega, sigmaomega

  !Read "ondes.txt" data  
  open (unit=15, file="ondes.txt", status='old')
  read(15,*)
  read(15,*)
  print*, '11'
  do i=1,42
    read (15,*) Aj,l(i),ls(i),F(i),D(i),Om(i),Me(i),Ve(i),Te(i),Ma(i),Ju(i), &
    Sa(i),Ur(i),Ne(i),Pa(i),Periode(i),ReREN(i),ImREN(i),ReMHB(i),ImMHB(i),ReCOR(i), &
    ImCOR(i),ReADD(i),ImADD(i),ReRET(i),ImRET(i)
        
    sigma(i) = l(i)*sigmal + ls(i)*sigmals + F(i)*sigmaF + D(i)*sigmaD + Om(i)*sigmaomega
    phi(i) = l(i)*phil + ls(i)*phils + F(i)*phiF + D(i)*phiD + Om(i)*phiomega
  end do
  print*,'12'
  close(15)

end program sum_error
