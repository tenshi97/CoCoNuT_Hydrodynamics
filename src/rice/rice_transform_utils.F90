module rice_transform_utils
  implicit none

contains

  subroutine rotate_vector_theta(u, u_new, sinbeta, cosbeta)
    real, intent(in)  :: u(0:3)
    real, intent(out) :: u_new(0:3)
    real, intent(in)  :: sinbeta, cosbeta

    u_new(0) = u(0)
    u_new(1) = cosbeta * u(1) - sinbeta * u(2)
    u_new(2) = cosbeta * u(2) + sinbeta * u(1)
    u_new(3) = u(3)

  end subroutine rotate_vector_theta

  subroutine rotate_vector_phi(u, u_new, sinbeta, cosbeta, sinpt, cospt)
    use rice_config, only: cartoon_grid
    real, intent(in)  :: u(0:3)
    real, intent(out) :: u_new(0:3)
    real, intent(in)  :: sinbeta, cosbeta
    real, intent(in)  :: sinpt, cospt

    if (.not. cartoon_grid) then
      u_new(0) = u(0)
      u_new(1) = cosbeta * u(1) + sinbeta * (-sinpt*u(3)) + cospt*u(1)*(1.0 - cosbeta)*cospt
      u_new(2) = cosbeta * u(2) - sinbeta * ( cospt*u(3)) + sinpt*u(2)*(1.0 - cosbeta)*sinpt
      u_new(3) = cosbeta * u(3) + sinbeta * ( sinpt*u(1)  + cospt*u(2))
    else
      u_new(0) = u(0)
      u_new(1) = cosbeta * u(1) - sinbeta * u(3)
      u_new(2) = cosbeta * u(2) + u(2)*(1.0 - cosbeta)
      u_new(3) = cosbeta * u(3) + sinbeta * u(1)
    endif

  end subroutine rotate_vector_phi

  subroutine rotate_vector_cartoon(u, u_new, sinbeta, cosbeta)
    real, intent(in)  :: u(0:3)
    real, intent(out) :: u_new(0:3)
    real, intent(in)  :: sinbeta, cosbeta

    u_new(0) = u(0)
    u_new(1) = u(1)
    u_new(2) = cosbeta * u(2) - sinbeta * u(3)
    u_new(3) = cosbeta * u(3) + sinbeta * u(2)

  end subroutine rotate_vector_cartoon

  subroutine get_nearest_index(target, x, n, isearch, period, ia, ib)
    real,    intent(in)  :: target
    integer, intent(in)  :: n
    real,    intent(in)  :: x(n)
    integer, intent(in)  :: isearch
    real,    intent(in)  :: period
    integer, intent(out) :: ia
    integer, intent(out) :: ib

    integer :: i

    logical :: found

    ! Returns i, the left bin (so i+1 is the right bin)
    ! Does a periodic search, which is not valid for mu,
    ! but mu should never cross boundary

    found = .false.

    if (isearch /= 1 .and. isearch /= n) then
      if (x(isearch-1) <= target .and. target <= x(isearch)) then
        ia = isearch-1
        ib = isearch
        found = .true.
      elseif (x(isearch) < target .and. target <= x(isearch+1)) then
        ia = isearch
        ib = isearch+1
        found = .true.
      endif
    elseif (isearch == 1) then ! left edge case
      if (x(1) <= target .and. target <= x(2)) then
        ia = 1
        ib = 2
        found = .true.
      elseif (x(n) <= target + period .and. target <= x(1)) then
        ia = n
        ib = 1
        found = .true.
      endif
    elseif (isearch == n) then ! right edge case
      if (x(n-1) <= target .and. target <= x(n)) then
        ia = n-1
        ib = n
        found = .true.
      elseif (x(n) <= target .and. target <= x(1) + period) then
        ia = n
        ib = 1
        found = .true.
      endif
    endif

    ! if quick search has failed (i.e. shifted by more than 1 cell)
    if (.not. found) then
      if (x(1) <= target .and. target <= x(n)) then
        do i = 1, n-1
          if (target <= x(i+1)) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      elseif (x(n) - period <= target .and. target <= x(1)) then
        ia = n
        ib = 1
      elseif (x(n) <= target .and. target <= x(1) + period) then
        ia = n
        ib = 1
      elseif (target <= x(n) - period) then
        do i = 1, n-1
          if (target <= x(i+1) - period) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      elseif (x(1) + period <= target) then
        do i = 1, n-1
          if (target <= x(i+1) + period) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      endif
    endif

  end subroutine get_nearest_index

  subroutine get_nearest_index_reverse(target, x, n, isearch, period, ia, ib)
    real,    intent(in)  :: target
    integer, intent(in)  :: n
    real,    intent(in)  :: x(n)
    integer, intent(in)  :: isearch
    real,    intent(in)  :: period
    integer, intent(out) :: ia
    integer, intent(out) :: ib

    integer :: i

    logical :: found

    ! Same algorith, but reversed signs

    found = .false.

    if (isearch /= 1 .and. isearch /= n) then
      if (-x(isearch-1) <= -target .and. -target <= -x(isearch)) then
        ia = isearch-1
        ib = isearch
        found = .true.
      elseif (-x(isearch) < -target .and. -target <= -x(isearch+1)) then
        ia = isearch
        ib = isearch+1
        found = .true.
      endif
    elseif (isearch == 1) then ! left edge case
      if (-x(1) <= -target .and. -target <= -x(2)) then
        ia = 1
        ib = 2
        found = .true.
      elseif (-x(n) <= -target + period .and. -target <= -x(1)) then
        ia = n
        ib = 1
        found = .true.
      endif
    elseif (isearch == n) then ! right edge case
      if (-x(n-1) <= -target .and. -target <= -x(n)) then
        ia = n-1
        ib = n
        found = .true.
      elseif (-x(n) <= -target .and. -target <= -x(1) + period) then
        ia = n
        ib = 1
        found = .true.
      endif
    endif

    ! if quick search has failed (i.e. shifted by more than 1 cell)
    if (.not. found) then
      if (-x(1) <= -target .and. -target <= -x(n)) then
        do i = 1, n-1
          if (-target <= -x(i+1)) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      elseif (-x(n) - period <= -target .and. -target <= -x(1)) then
        ia = n
        ib = 1
      elseif (-x(n) <= -target .and. -target <= -x(1) + period) then
        ia = n
        ib = 1
      elseif (-target <= -x(n) - period) then
        do i = 1, n-1
          if (-target <= -x(i+1) - period) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      elseif (-x(1) + period <= -target) then
        do i = 1, n-1
          if (-target <= -x(i+1) + period) then
            ia = i
            ib = i+1
            exit
          endif
        enddo
      endif
    endif

  end subroutine get_nearest_index_reverse

  function cartoon_rotation_angle(theta, alpha) result(angle)
    use rice_constants, only: pi

    real, intent(in) :: theta
    real, intent(in) :: alpha

    real :: cos_theta_new
    real :: e
    real :: angle

    cos_theta_new = cos(theta) * cos(alpha)
    e = asin(sin(theta) / sqrt(1.0 - cos_theta_new**2))
    angle = sign(0.5*pi - e, alpha * cos_theta_new)

  end function cartoon_rotation_angle

end module rice_transform_utils
