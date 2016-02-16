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

  complex(kind=xi),dimension(:),allocatable :: amplitude !amplitude of the terme at 457days 
  real(kind=xi), dimension(:), allocatable  :: ampl_time !corresponding time
  

  real(kind=xi) :: var !reading variable
  character(20) :: carc !reading variable

  real(kind=xi), dimension(2*Nbr_of_parameter) :: Ampl !complex amplitude that we need to adjust

  integer :: array_size !variable using to compute the size of the amplitude and ampl_time array
  integer :: mid_win, beg_win, end_win !variable to search the middle, the begin and the end of the sliding window
  integer :: beg_fit, end_fit !indice of the first and last possible midle of the window
  integer :: med_win !median indice of the sliding window

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
  !the window is sliding by a step to have the valeu of the term at many t(j) 
  !==============!


  !------------------------------!
  !allocation of the result array!
  !------------------------------!
  !array may be too big for the result but never too small (normaly)
  array_size = floor((t(Nbr_of_point)-t(1))/time_step)
  allocate(amplitude(array_size))
  allocate(ampl_time(array_size))
  ampl_time=0.0
  !------------------------------!

  !----------------------------------------------!
  !----------------------------------------------!
  !calculation of the amplitude at each time step!
  !----------------------------------------------!


  !looking for the first and last point where we can put the entire window around it
  i=1
  beg_fit_loop : do while (t(i)<(t(1)+Slide_Window/2))
    i=i+1
  end do beg_fit_loop
  beg_fit=i

  i=Nbr_of_point
  end_fit_loop : do while (t(i)>(t(Nbr_of_point)-Slide_Window/2))
    i=i-1
  end do end_fit_loop
  end_fit=i
  !---------------------------------------------------------------------------------
  ! print*,"array_size =", array_size
  ! print*,"beg_fit =",beg_fit
  ! print*,"end_fit =",end_fit

  j=1 !indice to travel in amplitude and ampl_time array
  mid_win=beg_fit
  amplitude_loop : do while (mid_win<=end_fit)
    i=mid_win

    end_loop : do while(t(i)<(t(mid_win)+Slide_Window/2)) ! looking for the end of sliding window
      i=i+1
      if (i==array_size) exit end_loop
    end do end_loop 
    end_win=i


    i=mid_win
    begin_loop : do while(t(i)>(t(mid_win) - Slide_Window/2 )) ! looking for the begening of the slidind window
      i=i-1
      if (i==1) exit begin_loop
    end do begin_loop
    beg_win=i


    !Porcessing the lsq in the previously calculate window
    call processing_lsq_period (period,Nbr_of_parameter,(end_win-beg_win+1),&
        dX(beg_win:end_win),dY(beg_win:end_win),errdX(beg_win:end_win),errdY(beg_win:end_win),&
        t(beg_win:end_win),&
        Ampl)
    amplitude(j)=cmplx(Ampl(1),Ampl(2))  
    !-----------------------------------------------------

    !looking for the median time of th sliding window
    med_win = (end_win - beg_win) / 2 + beg_win   
    ampl_time(j)=t(med_win)
    !------------------------------------------------
    ! print*,beg_win,mid_win,end_win,med_win,amplitude(j)

    i=mid_win+1 !avoid infinite loop
    midle_loop : do while(t(i)<t(mid_win)+time_step .and. i < Nbr_of_point) ! finding the next point 
      i=i+1
    end do midle_loop
    mid_win=i

    j=j+1
  end do amplitude_loop

  !----------------------------------------------!
  !----------------------------------------------!

  !------------!
  !writing data!
  !------------!
  open(11,file="457_days_ampl.dat",status="replace",&
    action="write",iostat=ios)
  
  do i = 1,array_size,1
    if(ampl_time(i)==0.0) exit
    write(11,"(3e26.16)") ampl_time(i), real(amplitude(i)), aimag(amplitude(i))
  end do

  close(11)
  !------------!


  deallocate(ampl_time)
  deallocate(amplitude)
  !==============!

end program lsq_457_days
