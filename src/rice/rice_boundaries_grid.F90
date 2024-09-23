module rice_boundaries_grid
  implicit none

contains

  subroutine fill_ghosts_grid(x, n, x_if, direction, cover_poles)
    integer,           intent(in)    :: n
    real,              intent(inout) :: x(-1:n+2)
    real,              intent(in)    :: x_if(-1:n+1)
    integer,           intent(in)    :: direction
    logical, optional, intent(in)    :: cover_poles

    logical :: cover_poles_switch

    if (present(cover_poles)) then
      cover_poles_switch = cover_poles
    else
      cover_poles_switch = .false.
    endif

    if (direction == 2 .and. cover_poles_switch) then ! interface area = 0 at poles
      ! this only affects reconstuction
      x(0)   = - x(2)
      x(-1)  = - x(3)
      x(n+1) = x_if(n) + (x_if(n) - x(n-1))
      x(n+2) = x_if(n) + (x_if(n) - x(n-2))
    elseif (direction == 2 .or. direction == 3) then ! periodic
      x(0)   = x_if(0) - (x_if(n) - x(n))
      x(-1)  = x(0) - (x(n) - x(n-1))
      x(n+1) = x_if(n) + (x(1) - x_if(0))
      x(n+2) = x(n+1) + (x(2) - x(1))
    elseif (direction == 1) then ! outflow
      x(0)   = x_if(0) - (x(1) - x_if(0))
      x(-1)  = x(0) - (x(2) - x(1)) ! mirrored
      x(n+1) = x_if(n) + (x_if(n) - x(n))
      x(n+2) = x(n+1) + (x(n+1) - x(n)) ! repeated
    endif
  end subroutine fill_ghosts_grid

  subroutine fill_ghosts_grid_if(x_if_g, n)
    integer, intent(in)    :: n
    real,    intent(inout) :: x_if_g(-1:n+1)

    x_if_g(-1) = x_if_g(0) - (x_if_g(1) - x_if_g(0))
    x_if_g(n+1) = x_if_g(n) + (x_if_g(n) - x_if_g(n-1))
  end subroutine fill_ghosts_grid_if

  subroutine fill_ghosts_scalar_1(x, n, direction)
    integer, intent(in)    :: n
    real,    intent(inout) :: x(0:n+1)
    integer, intent(in)    :: direction

    if (direction == 2 .or. direction == 3) then ! periodic
      x(0) = x(n)
      x(n+1) = x(1)
    elseif (direction == 1) then ! outflow
      x(0)  = x(1)
      x(n+1) = x(n)
    endif

  end subroutine fill_ghosts_scalar_1

end module rice_boundaries_grid
