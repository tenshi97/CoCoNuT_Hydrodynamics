module rice_boundaries
  use rice_grid, only: a_s, a_e, b_s, b_e, c_s, c_e

  implicit none

  interface fill_ghosts
   module procedure fill_ghosts_scalar, &
                    fill_ghosts_2d, &
                    fill_ghosts_3d, &
                    fill_ghosts_4d, &
                    fill_ghosts_f
  end interface

contains

  subroutine fill_boundary_ghosts
    use rice_config, only: nr, ntheta, nphi, neps, nmu, npsi, nflav
    use rice_grid,   only: f, f_eq, kappa_a, kappa_s, vfluid, alpha, beta, phiconf

    integer :: a, b, c

    !$omp parallel do
    do a = a_s, a_e
      do b = b_s, b_e
        call fill_ghosts(f_eq    (1:nflav,1:neps,                c_s-2:c_e+2,b,a), c_s, c_e, nphi, neps, nflav, 3)
        call fill_ghosts(kappa_a (1:nflav,1:neps,                c_s-2:c_e+2,b,a), c_s, c_e, nphi, neps, nflav, 3)
        call fill_ghosts(kappa_s (1:nflav,1:neps,                c_s-2:c_e+2,b,a), c_s, c_e, nphi, neps, nflav, 3)
        call fill_ghosts(vfluid  (1:3,                           c_s-2:c_e+2,b,a), c_s, c_e, nphi, 3,           3)
        call fill_ghosts(alpha   (                               c_s-2:c_e+2,b,a), c_s, c_e, nphi,              3)
        call fill_ghosts(beta    (1:3,                           c_s-2:c_e+2,b,a), c_s, c_e, nphi, 3,           3)
        call fill_ghosts(phiconf (                               c_s-2:c_e+2,b,a), c_s, c_e, nphi,              3)
        call fill_ghosts(f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c_s-2:c_e+2,b,a), c_s, c_e, nphi,              3, b)
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do a = a_s, a_e
      do c = c_s, c_e
        call fill_ghosts(f_eq    (1:nflav,1:neps,                c,b_s-2:b_e+2,a), b_s, b_e, ntheta, neps, nflav, 2)
        call fill_ghosts(kappa_a (1:nflav,1:neps,                c,b_s-2:b_e+2,a), b_s, b_e, ntheta, neps, nflav, 2)
        call fill_ghosts(kappa_s (1:nflav,1:neps,                c,b_s-2:b_e+2,a), b_s, b_e, ntheta, neps, nflav, 2)
        call fill_ghosts(vfluid  (1:3,                           c,b_s-2:b_e+2,a), b_s, b_e, ntheta, 3,           2)
        call fill_ghosts(alpha   (                               c,b_s-2:b_e+2,a), b_s, b_e, ntheta,              2)
        call fill_ghosts(beta    (1:3,                           c,b_s-2:b_e+2,a), b_s, b_e, ntheta, 3,           2)
        call fill_ghosts(phiconf (                               c,b_s-2:b_e+2,a), b_s, b_e, ntheta,              2)
        call fill_ghosts(f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b_s-2:b_e+2,a), b_s, b_e, ntheta,              2)
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do b = b_s, b_e
      do c = c_s, c_e
        call fill_ghosts(f_eq    (1:nflav,1:neps,                c,b,a_s-2:a_e+2), a_s, a_e, nr, neps, nflav, 1)
        call fill_ghosts(kappa_a (1:nflav,1:neps,                c,b,a_s-2:a_e+2), a_s, a_e, nr, neps, nflav, 1)
        call fill_ghosts(kappa_s (1:nflav,1:neps,                c,b,a_s-2:a_e+2), a_s, a_e, nr, neps, nflav, 1)
        call fill_ghosts(vfluid  (1:3,                           c,b,a_s-2:a_e+2), a_s, a_e, nr, 3,           1)
        call fill_ghosts(alpha   (                               c,b,a_s-2:a_e+2), a_s, a_e, nr,              1)
        call fill_ghosts(beta    (1:3,                           c,b,a_s-2:a_e+2), a_s, a_e, nr, 3,           1)
        call fill_ghosts(phiconf (                               c,b,a_s-2:a_e+2), a_s, a_e, nr,              1)
        call fill_ghosts(f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a_s-2:a_e+2), a_s, a_e, nr,              1)
      enddo
    enddo
    !$omp end parallel do

  end subroutine fill_boundary_ghosts

  subroutine fill_ghosts_scalar(x, i_s, i_e, n, direction)
    use rice_grid, only: grid_2d
    integer, intent(in)    :: i_s
    integer, intent(in)    :: i_e
    integer, intent(in)    :: n
    real,    intent(inout) :: x(i_s-2:i_e+2)
    integer, intent(in)    :: direction

    ! MPI has already handled periodic BCs

    if (direction == 2 .and. grid_2d) then ! repeated ghosts
      if (i_s-2 <=  -1) x( -1) = x(1)
      if (i_s-2 <=   0) x(  0) = x(1)
      if (i_e+2 >= n+1) x(n+1) = x(n)
      if (i_e+2 >= n+2) x(n+2) = x(n)
    elseif (direction == 1) then ! outflow
      if (i_s-2 <=  -1) x( -1) = x(1)
      if (i_s-2 <=   0) x(  0) = x(1)
      if (i_e+2 >= n+1) x(n+1) = x(n)
      if (i_e+2 >= n+2) x(n+2) = x(n)
    endif

  end subroutine fill_ghosts_scalar

  subroutine fill_ghosts_2d(x, i_s, i_e, n, n2, direction)
    use rice_grid, only: grid_2d
    integer, intent(in)    :: i_s
    integer, intent(in)    :: i_e
    integer, intent(in)    :: n
    integer, intent(in)    :: n2
    real,    intent(inout) :: x(1:n2,i_s-2:i_e+2)
    integer, intent(in)    :: direction

    ! MPI has already handled periodic BCs

    if (direction == 2 .and. grid_2d) then ! repeated ghosts
      if (i_s-2 <=  -1) x(:, -1) = x(:,1)
      if (i_s-2 <=   0) x(:,  0) = x(:,1)
      if (i_e+2 >= n+1) x(:,n+1) = x(:,n)
      if (i_e+2 >= n+2) x(:,n+2) = x(:,n)
    elseif (direction == 1) then ! outflow
      if (i_s-2 <=  -1) x(:, -1) = x(:,1)
      if (i_s-2 <=   0) x(:,  0) = x(:,1)
      if (i_e+2 >= n+1) x(:,n+1) = x(:,n)
      if (i_e+2 >= n+2) x(:,n+2) = x(:,n)
    endif

  end subroutine fill_ghosts_2d

  subroutine fill_ghosts_3d(x, i_s, i_e, n, n2, n3, direction)
    use rice_grid, only: grid_2d
    integer, intent(in)    :: i_s
    integer, intent(in)    :: i_e
    integer, intent(in)    :: n
    integer, intent(in)    :: n2
    integer, intent(in)    :: n3
    real,    intent(inout) :: x(1:n3,1:n2,i_s-2:i_e+2)
    integer, intent(in)    :: direction

    ! MPI has already handled periodic BCs

    if (direction == 2 .and. grid_2d) then ! repeated ghosts
      if (i_s-2 <=  -1) x(:,:, -1) = x(:,:,1)
      if (i_s-2 <=   0) x(:,:,  0) = x(:,:,1)
      if (i_e+2 >= n+1) x(:,:,n+1) = x(:,:,n)
      if (i_e+2 >= n+2) x(:,:,n+2) = x(:,:,n)
    elseif (direction == 1) then ! outflow
      if (i_s-2 <=  -1) x(:,:, -1) = x(:,:,1)
      if (i_s-2 <=   0) x(:,:,  0) = x(:,:,1)
      if (i_e+2 >= n+1) x(:,:,n+1) = x(:,:,n)
      if (i_e+2 >= n+2) x(:,:,n+2) = x(:,:,n)
    endif

  end subroutine fill_ghosts_3d

  subroutine fill_ghosts_4d(x, i_s, i_e, n, n2, n3, n4, direction)
    use rice_grid, only: grid_2d
    integer, intent(in)    :: i_s
    integer, intent(in)    :: i_e
    integer, intent(in)    :: n
    integer, intent(in)    :: n2
    integer, intent(in)    :: n3
    integer, intent(in)    :: n4
    real,    intent(inout) :: x(1:n4,1:n3,1:n2,i_s-2:i_e+2)
    integer, intent(in)    :: direction

    ! MPI has already handled periodic BCs

    if (direction == 2 .and. grid_2d) then ! repeated ghosts
      if (i_s-2 <=  -1) x(:,:,:, -1) = x(:,:,:,1)
      if (i_s-2 <=   0) x(:,:,:,  0) = x(:,:,:,1)
      if (i_e+2 >= n+1) x(:,:,:,n+1) = x(:,:,:,n)
      if (i_e+2 >= n+2) x(:,:,:,n+2) = x(:,:,:,n)
    elseif (direction == 1) then ! outflow
      if (i_s-2 <=  -1) x(:,:,:, -1) = x(:,:,:,1)
      if (i_s-2 <=   0) x(:,:,:,  0) = x(:,:,:,1)
      if (i_e+2 >= n+1) x(:,:,:,n+1) = x(:,:,:,n)
      if (i_e+2 >= n+2) x(:,:,:,n+2) = x(:,:,:,n)
    endif

  end subroutine fill_ghosts_4d

  subroutine fill_ghosts_f(x, i_s, i_e, n, direction, b)
    use rice_config,    only: neps, nmu, npsi, nflav, cartoon_grid, cartoon_rotate
    use rice_transform, only: transform_frame_cartoon
    use rice_Grid,      only: grid_2d

    integer,           intent(in)    :: i_s
    integer,           intent(in)    :: i_e
    integer,           intent(in)    :: n
    real,              intent(inout) :: x(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-2:i_e+2)
    integer,           intent(in)    :: direction
    integer, optional, intent(in)    :: b

    ! MPI has already handled periodic BCs

    if (direction == 3) then ! periodic
      if (cartoon_grid) then
        if (cartoon_rotate) then
          if (i_s == 1) x(:,:,:,:, -1) = transform_frame_cartoon(x(:,:,:,:,1), b, 1)
          if (i_s == 1) x(:,:,:,:,  0) = transform_frame_cartoon(x(:,:,:,:,1), b, 2)
          if (i_e == n) x(:,:,:,:,n+1) = transform_frame_cartoon(x(:,:,:,:,1), b, 3)
          if (i_e == n) x(:,:,:,:,n+2) = transform_frame_cartoon(x(:,:,:,:,1), b, 4)
        else
          if (i_s == 1) x(:,:,:,:, -1) = x(:,:,:,:,1)
          if (i_s == 1) x(:,:,:,:,  0) = x(:,:,:,:,1)
          if (i_e == n) x(:,:,:,:,n+1) = x(:,:,:,:,1)
          if (i_e == n) x(:,:,:,:,n+2) = x(:,:,:,:,1)
        endif
      endif

    elseif (direction == 2) then
      if (grid_2d) then ! repeated ghosts
        if (i_s-2 <=  -1) x(:,:,:,:, -1) = x(:,:,:,:,1)
        if (i_s-2 <=   0) x(:,:,:,:,  0) = x(:,:,:,:,1)
        if (i_e+2 >= n+1) x(:,:,:,:,n+1) = x(:,:,:,:,n)
        if (i_e+2 >= n+2) x(:,:,:,:,n+2) = x(:,:,:,:,n)
      endif

    elseif (direction == 1) then
      ! repeated
      if (i_s-2 <=  -1) x(:,:,:,:,-1) = x(:,:,:,:,1)
      if (i_s-2 <=   0) x(:,:,:,:, 0) = x(:,:,:,:,1)

      ! outflow (This is approximate, and does not take into account relativistic boosts)
      if (i_e+2 >= n+1) x(1:nflav,1:neps,1:npsi,1:nmu ,n+1) = x(1:nflav,1:neps,1:npsi,1:nmu,n)
      if (i_e+2 >= n+1) x(1:nflav,1:neps,1:npsi,-nmu:0,n+1) = 0.0
      if (i_e+2 >= n+2) x(1:nflav,1:neps,1:npsi,1:nmu ,n+2) = x(1:nflav,1:neps,1:npsi,1:nmu,n)
      if (i_e+2 >= n+2) x(1:nflav,1:neps,1:npsi,-nmu:0,n+2) = 0.0

    endif

  end subroutine fill_ghosts_f

end module rice_boundaries
