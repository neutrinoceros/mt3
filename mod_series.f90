Module mod_series

  implicit none
contains
  function ser(N,A_re,A_im,sigma,phi,time)

    implicit none

    integer, parameter        :: xi = selected_real_kind (15)

    complex(kind=xi) :: ser !return value of the function
    integer, intent(in) :: N !number of parameter
    real(kind=xi), dimension(N), intent(in) :: A_re, A_im !amplitude of the series
    real(kind=xi), dimension(N), intent(in) :: sigma, phi !frequence and phase of the amplitude
    real(kind=xi), intent(in) :: time !time to choose for the calculation
    real(kind=xi) :: X, Y

    integer :: j

    X = 0._xi
    Y = 0._xi

    !the 
    do j=1,N
      X = A_re(j)*cos(sigma(j)*time+phi(j)) - A_im(j)*sin(sigma(j)*time+phi(j))    
      Y = A_re(j)*sin(sigma(j)*time+phi(j)) + A_im(j)*cos(sigma(j)*time+phi(j))    
    end do

    ser = cmplx(X, Y)

  end function ser

end module mod_series
