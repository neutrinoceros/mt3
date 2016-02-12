program lsq_457_days

  !program in order to fit the complex amplitude of the term at 450 day

  use arg_nut
  use mod_matrix
  use mod_lsq
  implicit none


  integer,parameter        :: Nbr_of_point=5980,Nbr_of_parameter=1
  real (kind=xi),parameter :: period=457._xi!days period wich we want to fit the complex amplitude
  real (kind=xi),parameter :: time_step =365.
  real (kind=xi),parameter :: Slide_Window = 2.*365. !days, it is the windows where we process the amplitude

  real(kind=xi), dimension(Nbr_of_point)   :: t, dX, dY, errdX, errdY, corrdXdY !parameter of observed nutation

  real(kind=xi), dimension(:), allocatable :: amplitude, ampl_time !amplitude of the terme at 457days and corresponding time
  real(kind=xi) :: date
  

  real(kind=xi) :: var !reading variable
  character(20) :: carc !reading variable

  real(kind=xi), dimension(2*Nbr_of_parameter) :: Ampl !complex amplitude that we need to adjust

  integer :: array_size !variable using to compute the size of the amplitude and ampl_time array
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

  !==============!
  !sliding window!
  !--------------!
  !we use a sliding window to fit the amplitude of 457 days terme.
  !the window is sliding by a step to have the valeu of the term at many date 
  !==============!


  !------------------------------!
  !allocation of the result array!
  !------------------------------!
  !array may be too big for the result but never too small (normaly)
  array_size = floor((t(Nbr_of_point)-t(1))/time_step)
  allocate(amplitude(array_size))
  allocate(ampl_time(array_size))
  !------------------------------!

  !==============================================!
  !calculation of the amplitude at each time step!
  !----------------------------------------------!

  amplitude_loop : do j=1,array_size,1
    i=j
    date=t(j)+time_step
    do while(t(i)<date) ! finding the next point 
      i=i+1

      if (i>array_size) exit amplitude_loop
    end do 

  end do amplitude_loop



  !==============================================!

  ! call processing_lsq_period(period,Nbr_of_parameter,Nbr_of_point,&
  !    dX,dY,errdX,errdY,t,&
  !    Ampl)

  deallocate(ampl_time)
  deallocate(amplitude)
end program lsq_457_days
