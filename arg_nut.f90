module arg_nut

!======================================================================    
! script used to define the Fundamental Argument of Nutation (in radian)
!======================================================================    

  implicit none

!======================================================================  
! Variables declaration
!======================================================================  
  
  integer, parameter:: xi = selected_real_kind (15)
  real (kind=xi), parameter :: pi=4._xi*atan(1._xi)

  !phi.. (phase) equal to constant part of the equation
  real(kind=xi), parameter :: phil = 134.96340251_xi*pi/180._xi 
  real(kind=xi), parameter :: phils = (357.52910918_xi*pi/180._xi)
  real(kind=xi), parameter :: phiF = 93.27209062_xi*pi/180._xi
  real(kind=xi), parameter :: phiD = (297.85019547_xi*pi/180._xi)
  real(kind=xi), parameter :: phiomega = 125.04455501_xi*pi/180._xi

  !sigma.. (frequency) equal to coeficient of "t" of the equation
  real(kind=xi), parameter :: sigmal =(1717915923.2178_xi/3600._xi)*pi/180._xi
  real(kind=xi), parameter :: sigmals =(129596581.0481_xi/3600._xi)*pi/180._xi
  real(kind=xi), parameter :: sigmaF =(1739527262.8478_xi/3600._xi)*pi/180._xi
  real(kind=xi), parameter :: sigmaD =(1602961601.2090_xi/3600._xi)*pi/180._xi
  real(kind=xi), parameter :: sigmaomega =(-6962890.5431_xi/3600._xi)*pi/180._xi

  private 

  public :: xi, pi, phil, phils, phiF, phiD, phiomega, sigmal, sigmals, sigmaF, sigmaD, sigmaomega

contains
end module arg_nut
