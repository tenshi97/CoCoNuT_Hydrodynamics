 module rice_laxwendroff
  implicit none

contains

  subroutine lw_state(f_lw, f_l_if, f_r_if, &
                      x_if, x_l, x_r, &
                      f_eq_l, f_eq_r, kappa_a_l, kappa_a_r, kappa_s_l, kappa_s_r, &
                      volume_l, volume_r, area, &
                      dt, diffusive_limit, direction)
    use rice_config,     only: neps, nmu, npsi, nflav
    use rice_phase,      only: phase_step
    use rice_grid,       only: u_com
    use rice_flux,       only: flux_mono

    real,     intent(out) :: f_lw  (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: f_l_if(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: f_r_if(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: x_if
    real,     intent(in)  :: x_l
    real,     intent(in)  :: x_r
    real,     intent(in)  :: f_eq_l   (1:nflav,1:neps)
    real,     intent(in)  :: f_eq_r   (1:nflav,1:neps)
    real,     intent(in)  :: kappa_a_l(1:nflav,1:neps)
    real,     intent(in)  :: kappa_a_r(1:nflav,1:neps)
    real,     intent(in)  :: kappa_s_l(1:nflav,1:neps)
    real,     intent(in)  :: kappa_s_r(1:nflav,1:neps)
    real,     intent(in)  :: volume_l
    real,     intent(in)  :: volume_r
    real,     intent(in)  :: area
    real,     intent(in)  :: dt
    logical,  intent(in)  :: diffusive_limit
    integer,  intent(in)  :: direction

    real :: flux_left (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_right(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_if   (1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: u_com_d(1:npsi,-nmu:nmu)

    real :: f_eq_if   (1:nflav,1:neps)
    real :: kappa_a_if(1:nflav,1:neps)
    real :: kappa_s_if(1:nflav,1:neps)
    real :: meanvol

    real :: w_l, w_r

    call interpolation_weight(x_l, x_r, x_if, w_l, w_r)

    ! Interpolate between the two neighbouring cells
    f_lw = f_l_if*w_l + f_r_if*w_r

    ! Calculate flux from left and right cells to advance intermediate state
    ! Make an approximation that these cells are in the same frame (use ucom in flux call)

    u_com_d = u_com(direction,:,:)

    ! Flux from the left cell
    flux_left = flux_mono(f_l_if, u_com_d)
    ! Flux from the right cell
    flux_right = flux_mono(f_r_if, u_com_d)

    ! Interpolate source terms
    f_eq_if    = f_eq_l   *w_l + f_eq_r   *w_r
    kappa_a_if = kappa_a_l*w_l + kappa_a_r*w_r
    kappa_s_if = kappa_s_l*w_l + kappa_s_r*w_r

    ! Average volume between the two cells
    meanvol = volume_l*w_l + volume_r*w_r

    ! Add fluxes together (using the interface area, and average volume)
    ! These can add because we approximate them to be in the same frame
    flux_if = (flux_left - flux_right) * area / meanvol

    ! Do the implicit solve to get the half time step
    call phase_step(f_lw, flux_if, &
                    f_eq_if   (1:nflav,1:neps), &
                    kappa_a_if(1:nflav,1:neps), &
                    kappa_s_if(1:nflav,1:neps), &
                    0.5*dt, &
                    diffusive_limit)

  end subroutine lw_state

  subroutine interpolation_weight(x_l, x_r, x_if, w_l, w_r)
    real, intent(in)  :: x_l
    real, intent(in)  :: x_r
    real, intent(in)  :: x_if
    real, intent(out) :: w_l
    real, intent(out) :: w_r

    w_r = (x_if - x_l) / (x_r - x_l)
    w_l = 1.0 - w_r
  end subroutine interpolation_weight

end module rice_laxwendroff
