module rice_dt
  implicit none

contains

  subroutine get_dt(dt)
    use rice_config, only: nr, ntheta, nphi, cfl_factor, clight, average_poles
    use rice_grid, only: volume, area_phi, area_theta, area_r

    real, intent(out) :: dt

    real :: dx
    real :: dxmin
    real :: dr
    real :: dtheta
    real :: dphi

    integer :: b0, b1

    integer :: a, b, c
    integer :: limiting_direction

    dxmin = 1.0e99

    ! Allow larger time step when averaging poles
    if (average_poles) then
      b0 = 2
      b1 = ntheta-1
    else
      b0 = 1
      b1 = ntheta
    endif

    do a = 1, nr
      do b = b0, b1
        do c = 1, nphi
          dr      = volume(c,b,a) / max(area_r(c,b,a), area_r(c,b,a-1))
          dtheta  = volume(c,b,a) / max(area_theta(c,b,a), area_theta(c,b-1,a))
          dphi    = volume(c,b,a) / max(area_phi(c,b,a), area_phi(c-1,b,a))

          dx = min(dr, min(dtheta, dphi))
          if (dx < dxmin) then
            dxmin = dx
            if (dr <= dtheta .and. dr <= dphi) then
              limiting_direction = 1
            elseif (dtheta <= dphi) then
              limiting_direction = 2
            else
              limiting_direction = 3
            endif
          endif
        enddo
      enddo
    enddo

    print*, '(dt) Limiting direction:', limiting_direction, 'dxmin =', dxmin

    dt = cfl_factor * dxmin / clight

    print*, '(dt) dt =', dt

  end subroutine get_dt

end module rice_dt
