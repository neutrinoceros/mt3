program lsq_450_days

  !program in order to fit the complex amplitude of the term at 450 day

  use arg_nut
  use mod_matrix
  implicit none


  integer,parameter        :: Nbr_of_point=5980,Nbr_of_parameter=1
  real (kind=xi),parameter :: period=450._xi!days period wich we want to fit the complex amplitude
  real (kind=xi)           :: sigma = 2._xi*pi/period

  real(kind=xi), dimension(Nbr_of_point)   :: t, dX, dY, errdX, errdY, corrdXdY !parameter of observed nutation
  real(kind=xi), dimension(2*Nbr_of_point) :: dXdY !matrix of observable


  real(kind=xi) :: var !reading variable
  character(20) :: carc !reading variable

  real(kind=xi), dimension(2*Nbr_of_point,2*Nbr_of_parameter) :: M !matrix like dXdY = M Ampl
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_parameter) :: MM, P !MM is equal to tranpose(M) matrix product M
  real(kind=xi), dimension(2*Nbr_of_parameter,2*Nbr_of_point) :: Q, MMM
  real(kind=xi), dimension(2*Nbr_of_parameter) :: Ampl !complex amplitude that we need to adjust
  real(kind=xi), dimension(Nbr_of_parameter) :: A, B
  integer :: i, j, k, r, s !loop variable
  integer :: ios !checking io variable

  !Read "opa2015a.eops" data  
  open (unit=10, file="data/opa2015a.eops", status='old',&
    iostat=ios)
  if (ios/=0) stop "pb ouverture opa2015a.eops"

  do k=1,38!skipping the header
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
       errdX(j) = 100000.
    end if
    if (errdY(j) .lt. 0.0009) then
       errdY(j) = 100000.
    end if

  end do
  close(10)

  !Creation of the matrix M
  do j=1,Nbr_of_point
    !processing the 1st cadran of the matrix
    ! do k=1,Nbr_of_parameter
    !   M(j,k) = 1./errdX(j)*cos(sigma(k)*t(j)+phi(k))
    ! end do

    !processing the 2nd cadran of the matrix
    ! do k=1,Nbr_of_parameter
    !   r = Nbr_of_parameter + k
    !   M(j,r) = -1./errdX(j)*sin(sigma(k)*t(j)+phi(k))
    ! end do
    M(j,1) = 1./errdX(j)*cos(sigma*t(j))
    M(j,2) = 1./errdX(j)*sin(sigma*t(j))

  end do

  do j = 1,Nbr_of_point
    s = Nbr_of_point + j

    !processing the 3rd cadran of the matrix
    ! do k=1,Nbr_of_parameter
    !   M(s,k) = 1./errdY(j)*sin(sigma(k)*t(j)+phi(k))
    ! end do

    !processing the 4th cadran of the matrix
    ! do k=1,Nbr_of_parameter
    !   r = Nbr_of_parameter + k
    !   M(s,r) = 1./errdY(j)*cos(sigma(k)*t(j)+phi(k))
    ! end do
    M(s,1) = -1./errdY(j)*sin(sigma*t(j))
    M(s,2) = 1./errdY(j)*cos(sigma*t(j))

  end do




  
 
end program lsq_450_days
