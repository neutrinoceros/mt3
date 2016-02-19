module mod_lsq
  use mod_matrix
  use arg_nut
  implicit none
contains

  !module to compute 
  subroutine processing_lsq_period (periode,Nbr_of_parameter,Nbr_of_point,&
      measure_array_X,measure_array_Y,error_array_X,error_array_Y,time,&
      amplitude_array,error_amplitude)
    !Subroutine to processe the amplitude of 457 days nutation periode
    implicit none

    !========================!
    !argument of the function!
    !------------------------!
    real(kind=xi),                         intent (in) :: periode ! periode that we want to fit the parameter
    integer,                               intent (in) :: Nbr_of_parameter ! nbr of parameter we want to fit
    integer,                               intent (in) :: Nbr_of_point ! nbr of measure
    real(kind=xi),dimension(Nbr_of_point), intent (in) :: measure_array_X, measure_array_Y ! array containing the measure point
    real(kind=xi),dimension(Nbr_of_point), intent (in) :: error_array_X,error_array_Y !array contains the error of measure
    real(kind=xi),dimension(Nbr_of_point), intent (in) :: time !time of th measure array

    real(kind=xi),dimension(2*Nbr_of_parameter), intent(out) :: amplitude_array ! amplitude we want to process 
    real(kind=xi),dimension(2*Nbr_of_parameter), intent(out) :: error_amplitude !error in the amplitude processing by lsq
    !========================!

    !==============!
    !local variable!
    !--------------!

    real(kind=xi) :: sigma 

    real(kind=xi), dimension(2*Nbr_of_point) :: dXdY !matrix of observable

    real(kind=xi), dimension(2*Nbr_of_point,2*Nbr_of_parameter)     :: M !matrix like dXdY = M Ampl
    real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_parameter) :: MM, P !MM is equal to tranpose(M) matrix product M
    real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_point)     :: Q, MMM
    real(kind=xi), dimension(2*Nbr_of_parameter)                    :: Ampl !complex amplitude that we need to adjust
    real(kind=xi), dimension(2*Nbr_of_point) :: MA, tmpMa !matrix to compute the error
    real(kind=xi) ::sigma_lsq ! residu

    integer :: i, j, k, r, s !loop variable

    !==============!
    sigma = 2._xi*pi/periode

    !Creation of the matrix M
    do j=1,Nbr_of_point
      s = Nbr_of_point + j

      M(j,1) = 1./error_array_X(j)*cos(sigma*time(j))
      M(j,2) = 1./error_array_X(j)*sin(sigma*time(j))

      M(s,1) = -1./error_array_Y(j)*sin(sigma*time(j))
      M(s,2) = 1./error_array_Y(j)*cos(sigma*time(j))

    end do

    !computation of the observationnal array
    dXdY(:Nbr_of_point)   = measure_array_X/error_array_x
    dXdY(Nbr_of_point+1:) = measure_array_Y/error_array_Y

    !Inverse matrix with A = ((M[t]*M)^-1)*(M[t]*X)
    Q = transpose(M) ! Transopse the Matrix
    MM = matmul(Q,M) ! Multiply the Matrix
    P = inv(MM)      ! Inverse the Matrix using function which defined
    MMM = matmul(P,Q)

    amplitude_array = matmul(MMM,dXdY)

    !processing the error
    MA    = matmul(M,amplitude_array)
    tmpMa = (dXdY-MA)**2
    sigma_lsq = sqrt(sum(tmpMa))/(Nbr_of_point-Nbr_of_parameter)
    
    do i = 1, 2*Nbr_of_parameter,1
      error_amplitude(i) = sqrt(P(i,i))
      ! print*, error_amplitude(i)
    end do

    error_amplitude=error_amplitude*sigma_lsq


  end subroutine processing_lsq_period

end module mod_lsq

