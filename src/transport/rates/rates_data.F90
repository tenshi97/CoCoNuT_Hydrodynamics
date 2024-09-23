module rates_data_type

! New module for handling all interaction rates that are tabulated
!
! A. Marek, Jul. 2009, MPA

  use precision

  implicit none

  public lms_rates_table

  save

  type rates_tables
     ! this data type is used for various rates such as
  ! the LMS-electron capture rates and the ISN-rates
     ! feel free to add more

     ! variables that are used in all tables

     character(len=30) :: name

     integer(kind=ik) :: nro, ntt, nye, nene
     real(kind=rk)    :: romax, romin, ttmax, ttmin, yemax, yemin


     real(kind=rk), dimension(:),   pointer   :: rho, tem,log_rho, energy
     real(kind=rk), dimension(:,:), pointer   :: ye

     real(kind=rk), dimension(:,:,:), pointer :: chel, yp, chep, yn, chen, &
                                                 zav, aav, zh, ah, yhvy,   &
                                                 yalpha, rtot, rnrm, qeff, &
                                                 mu_nu, mu_e, n

     real(kind=rk), dimension(:,:,:,:), pointer :: spectra


     real(kind=rk), dimension(:,:,:,:,:), pointer :: scatter


  end type rates_tables

  type(rates_tables) :: lms_rates_table(1)

  type(rates_tables) :: isn_rates_table(1)

end module rates_data_type
