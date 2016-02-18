Module mod_series
  use arg_nut
  implicit none
contains


  complex (kind=xi) function ser(N,A_re,A_im,sigma,phi,time)
    !Function to compute the value of the nutation with the complex coeficient A (which is pass by A_re and A_im)
    !and sigma and phi at the time t

    implicit none

    !------------------------!
    !argument of the function!
    !------------------------!
    integer,                     intent (in) :: N          ! number of parameter
    real(kind=xi), dimension(N), intent (in) :: A_re, A_im ! amplitude of the series
    real(kind=xi), dimension(N), intent (in) :: sigma, phi ! frequence and phase of the amplitude
    real(kind=xi),               intent (in) :: time       ! time to choose for the calculation
    !------------------------!

    !--------------!
    !local variable!
    !--------------!
    real(kind=xi) :: X, Y ! temporar variable in order to compute the real and the imaginary part of the nutation
    integer       :: j    ! loop variable
    !--------------!

    X = 0._xi
    Y = 0._xi

    !the 
    do j=1,N
      X = X + A_re(j)*cos(sigma(j)*time+phi(j)) - A_im(j)*sin(sigma(j)*time+phi(j))    
      Y = Y + A_re(j)*sin(sigma(j)*time+phi(j)) + A_im(j)*cos(sigma(j)*time+phi(j))    
    end do

    ser = cmplx(X, Y)

  end function ser

end module mod_series
