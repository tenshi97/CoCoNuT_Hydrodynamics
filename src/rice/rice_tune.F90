module rice_tune
  implicit none

contains

subroutine tune_areas(preset)
  use rice_config, only: nr, ntheta, nphi, lw_only
  use rice_grid,   only: f, kappa_a, kappa_s, f_eq, vfluid, phiconf, alpha, beta, &
                         area_theta, area_phi
  use rice_dt,     only: get_dt
  use rice_step,   only: step

  real, optional, intent(in) :: preset

  real :: dt

  real :: xtune(nr)
  real :: factor(nr)
  real :: fh_old(nr)
  real :: fh_new
  real :: maxerr
  real :: area_theta_org(1:nphi,0:ntheta,1:nr)
  real :: area_phi_org  (0:nphi,1:ntheta,1:nr)

  integer :: a, iit

  print*, '(tune) Beginning area tune-up'

  area_theta_org = area_theta
  area_phi_org   = area_phi

  if (.not. present(preset)) then
    call get_dt(dt)

    print*, '(tune) Setting tuning ICs'
    kappa_a = 0.0
    kappa_s = 0.0
    f_eq = 0.0
    vfluid = 0.0
    phiconf = 1.0
    alpha = 1.0
    beta = 0.0

    xtune = 1.0
    factor = 50.0
    fh_old = 1.0

    lw_only = .true.

    print*, '(tune) Performing iterations'
    do iit = 1, 10000

      f = 1.0
      call step(dt)

      maxerr = 0.0
      do a = 1, nr-1

        fh_new = fh(f(:,:,:,:,1,1,a))
        maxerr = max(maxerr, abs(fh_new))

        if (fh_new / fh_old(a) < 0.0) then
          factor(a) = factor(a) * 0.5
        endif

        fh_old(a) = fh_new

        xtune(a) = xtune(a) - factor(a) * fh_new

        area_theta(1,0,a) = area_theta_org(1,0,a) * xtune(a)
        area_theta(1,1,a) = area_theta_org(1,1,a) * xtune(a)
        area_phi  (0,1,a) = area_phi_org  (0,1,a) * xtune(a)
        area_phi  (1,1,a) = area_phi_org  (1,1,a) * xtune(a)
      enddo

      ! print*, iit, maxerr
      if (maxerr < 1.e-16) exit

    enddo

    ! do a = 1, nr
    !   print*, a, xtune(a), factor(a), fh_old(a)
    ! enddo

    lw_only = .false.

  else
    print*, '(tune) Using tuning preset =', preset
    xtune(:) = preset
  endif

  do a = 1, nr
    area_theta(1,0,a) = area_theta_org(1,0,a) * xtune(a)
    area_theta(1,1,a) = area_theta_org(1,1,a) * xtune(a)
    area_phi  (0,1,a) = area_phi_org  (0,1,a) * xtune(a)
    area_phi  (1,1,a) = area_phi_org  (1,1,a) * xtune(a)
  enddo

  print*, '(tune) Tune-up complete, iteration number =', iit
  print*, '(tune) Max H error =', maxerr

end subroutine tune_areas

function fj(f) result(j_out)
  use rice_constants, only: pi
  use rice_config, only: nmu, npsi, neps, nflav
  use rice_grid, only: domega
  real, intent(in) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)

  real :: j_out
  integer :: j, k

  j_out = 0.0
  do j = -nmu, nmu
    do k = 1, npsi
      j_out = j_out + f(1,1,k,j) * domega(k,j)
    enddo
  enddo

  j_out = j_out / (4.0 * pi)

end function fj

function fh(f) result(h_out)
  use rice_constants, only: pi
  use rice_config, only: nmu, npsi, neps, nflav
  use rice_grid, only: u_com, domega
  real, intent(in) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)

  real :: h_out
  integer :: j, k

  h_out = 0.0
  do j = -nmu, nmu
    do k = 1, npsi
      h_out = h_out + f(1,1,k,j) * u_com(1,k,j) * domega(k,j)
    enddo
  enddo

  h_out = h_out / (4.0 * pi)

end function fh


end module rice_tune
