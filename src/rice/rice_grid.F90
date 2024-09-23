module rice_grid
  use rice_config, only: nr, ntheta, nphi, neps, nmu, npsi, nflav, heun
  implicit none

  ! Grid interfaces
  real, allocatable :: r_if     (:)
  real, allocatable :: theta_if (:)
  real, allocatable :: phi_if   (:)
  real, allocatable :: eps_if   (:)
  real, allocatable :: mu_if    (:)
  real, allocatable :: psi_if   (:)

  ! Grid centers
  real, allocatable :: r     (:)
  real, allocatable :: theta (:)
  real, allocatable :: phi   (:)
  real, allocatable :: eps   (:)
  real, allocatable :: mu    (:)
  real, allocatable :: psi   (:)

  ! Jump switch
  logical, allocatable :: jump (:,:,:)

  ! Auxilliary variables
  real, allocatable :: cpsi(:)
  real, allocatable :: spsi(:)

  ! Grid cell volume
  real, allocatable :: volume(:,:,:)

  ! Grid cell areas
  real, allocatable :: area_r     (:,:,:)
  real, allocatable :: area_theta (:,:,:)
  real, allocatable :: area_phi   (:,:,:)

  ! Equivalent volumes and areas in spherical coordinates (when using cartoon grid)
  real, allocatable :: volume_spherical    (:,:,:)
  real, allocatable :: area_r_spherical    (:,:,:)
  real, allocatable :: area_theta_spherical(:,:,:)
  real, allocatable :: area_phi_spherical  (:,:,:)

  ! Fraction of radial area contribution from the inner radial interface
  real, allocatable :: cell_ratio (:,:,:)

  ! Phase space area
  real :: omega
  real, allocatable :: domega        (:,:)
  real, allocatable :: domega_scaled (:,:)
  real, allocatable :: domega_inv    (:,:)

  ! Phase space volume
  real, allocatable :: epsvol     (:)
  real, allocatable :: epsvol_inv (:)

  ! Differentials
  real, allocatable :: eps2domegadeps   (:,:,:)
  real, allocatable :: eps3domegadeps   (:,:,:)
  real, allocatable :: mueps3domegadeps (:,:,:)
  real, allocatable :: mu2eps3domegadeps(:,:,:)

  ! Radiation velocity
  real, allocatable :: u_com(:,:,:)       ! same everywhere
  real, allocatable :: u_lab(:,:,:,:,:,:) ! depends on vfluid, alpha, beta, and phiconf

  ! Fluid velocity
  real, allocatable :: vfluid    (:,:,:,:)
  real, allocatable :: vfluid_old(:,:,:,:)

  ! f arrays
  real, allocatable :: f     (:,:,:,:,:,:,:)
  real, allocatable :: f_old (:,:,:,:,:,:,:)
  real, allocatable :: f_null(:,:,:,:      ) ! a single spatial cell of zeroes

  ! Source term arrays
  real, allocatable :: f_eq   (:,:,:,:,:)
  real, allocatable :: kappa_a(:,:,:,:,:)
  real, allocatable :: kappa_s(:,:,:,:,:)

  ! Time
  real :: time = 0.0

  ! Output index
  integer :: index = 0

  ! Number of steps
  integer :: nstep = 0

  ! GR terms
  real, allocatable :: alpha      (:,:,:) ! lapse
  real, allocatable :: alpha_old  (:,:,:)
  real, allocatable :: beta     (:,:,:,:) ! shift
  real, allocatable :: beta_old (:,:,:,:)
  real, allocatable :: phiconf    (:,:,:) ! conformal factor
  real, allocatable :: phiconf_old(:,:,:)

  ! Cartoon rotation weights
  real,    allocatable :: cartoon_weight (:,:,:)
  integer, allocatable :: cartoon_index(:,:,:,:)

  ! Communication buffer
  real, target, allocatable :: sbuf(:) ! send
#ifdef MPI_BT
  real, allocatable :: rbuf(:) ! receive
#else
  real, pointer :: rbuf(:)
#endif /* MPI_BT */
  integer :: nbuf       ! max buffer size
  integer :: ibuf       ! current buffer size
  integer :: ibuf_send  ! size of buffer sent

  ! Indicator for 2d runs
  logical :: grid_2d = .false.

  ! Switch to force externally provided coordinates to be used, and prevent optimisations
  logical :: external_coordinates = .false.

  ! Domain indices
  !          r   theta phi
  integer :: a_s, b_s, c_s ! start
  integer :: a_e, b_e, c_e ! end
  integer :: na,  nb,  nc  ! size on this task

  ! MPI coordinates
  integer :: cart_coords(3)

contains

  subroutine setup_grid

    call set_domain
    call allocate_memory
    call spatial_grids
    call phase_grids
    call aux_grids

  end subroutine setup_grid

  subroutine set_domain
    ! Choose the part of the grid that this thread solves
    ! For non-MPI runs, set the domain to the entire grid
#ifdef MPI_BT
    use rice_config, only: nr, ntheta, nphi, mpi_nr, mpi_ntheta, mpi_nphi

    if (mod(nr, mpi_nr) /= 0) stop 'nr not divisible by mpi_nr'
    if (mod(ntheta, mpi_ntheta) /= 0) stop 'ntheta not divisible by mpi_ntheta'
    if (mod(nphi, mpi_nphi) /= 0) stop 'nphi not divisible by mpi_nphi'

    na = nr / mpi_nr
    a_s = cart_coords(1)*na + 1
    a_e = cart_coords(1)*na + na

    nb = ntheta / mpi_ntheta
    b_s = cart_coords(2)*nb + 1
    b_e = cart_coords(2)*nb + nb

    nc = nphi / mpi_nphi
    c_s = cart_coords(3)*nc + 1
    c_e = cart_coords(3)*nc + nc

#else

    a_s = 1
    a_e = nr
    na = a_e - a_s + 1

    b_s = 1
    b_e = ntheta
    nb = b_e - b_s + 1

    c_s = 1
    c_e = nphi
    nc = c_e - c_s + 1
#endif

  end subroutine set_domain

  subroutine allocate_memory

    integer :: allocstat

    ! -- Constants throughout the entire run --

    !  g | g | . | ... | . | g | g
    !   -1   0   1     n-1 n   n+1    -- interface indexing
    ! -1   0   1         n  n+1  n+2  -- center indexing
    !
    ! . : cell center
    ! | : interface
    ! g : ghost center
    !
    ! Momentum coordinates do not need ghost zones

    ! Interface coordinates
    allocate(r_if     (-1:nr+1),      stat=allocstat)
    allocate(theta_if (-1:ntheta+1),  stat=allocstat)
    allocate(phi_if   (-1:nphi+1),    stat=allocstat)
    allocate(eps_if   (0:neps),       stat=allocstat)
    allocate(mu_if    (-nmu-1:nmu),   stat=allocstat)
    allocate(psi_if   (0:npsi),       stat=allocstat)
    r_if     = 0.0
    theta_if = 0.0
    phi_if   = 0.0
    eps_if   = 0.0
    mu_if    = 0.0
    psi_if   = 0.0

    ! Cell center coordinates
    allocate(r     (-1:nr+2),         stat=allocstat)
    allocate(theta (-1:ntheta+2),     stat=allocstat)
    allocate(phi   (-1:nphi+2),       stat=allocstat)
    allocate(eps   (1:neps),          stat=allocstat)
    allocate(mu    (-nmu:nmu),        stat=allocstat)
    allocate(psi   (1:npsi),          stat=allocstat)
    r     = 0.0
    theta = 0.0
    phi   = 0.0
    eps   = 0.0
    mu    = 0.0
    psi   = 0.0

    ! Pre-computed trig values (constant)
    allocate(cpsi(1:npsi),        stat=allocstat)
    allocate(spsi(1:npsi),        stat=allocstat)
    cpsi = 0.0
    spsi = 0.0

    ! Cell volumes
    allocate(volume(0:nphi+1,0:ntheta+1,0:nr+1), stat=allocstat)
    volume = 0.0

    ! Cell surface areas
    allocate(area_r(1:nphi,1:ntheta,0:nr),       stat=allocstat)
    allocate(area_theta(1:nphi,0:ntheta,1:nr),   stat=allocstat)
    allocate(area_phi(0:nphi,1:ntheta,1:nr),     stat=allocstat)
    area_r     = 0.0
    area_theta = 0.0
    area_phi   = 0.0

    ! Surface element of the unit sphere for momentum space
    allocate(domega(1:npsi,-nmu:nmu),            stat=allocstat)
    allocate(domega_scaled(1:npsi,-nmu:nmu),     stat=allocstat)
    allocate(domega_inv(1:npsi,-nmu:nmu),        stat=allocstat)
    domega        = 0.0
    domega_scaled = 0.0
    domega_inv    = 0.0

    ! Volume in energy space
    allocate(epsvol(1:neps),                     stat=allocstat)
    allocate(epsvol_inv(1:neps),                 stat=allocstat)
    epsvol     = 0.0
    epsvol_inv = 0.0

    ! Differentials for the integrals
    allocate(eps2domegadeps   (1:neps,1:npsi,-nmu:nmu), stat=allocstat)
    allocate(eps3domegadeps   (1:neps,1:npsi,-nmu:nmu), stat=allocstat)
    allocate(mueps3domegadeps (1:neps,1:npsi,-nmu:nmu), stat=allocstat)
    allocate(mu2eps3domegadeps(1:neps,1:npsi,-nmu:nmu), stat=allocstat)
    eps2domegadeps    = 0.0
    eps3domegadeps    = 0.0
    mueps3domegadeps  = 0.0
    mu2eps3domegadeps = 0.0

    ! Comoving four-velocities
    allocate(u_com(0:3,1:npsi,-nmu:nmu),        stat=allocstat)
    u_com = 0.0

    ! Ratio of radial flux travelling through angular interfaces (diagnostic)
    allocate(cell_ratio(1:nphi,1:ntheta,1:nr),    stat=allocstat)
    cell_ratio = 0.0


    ! -- Variables --

    ! Lab-frame four-velocities
    allocate(u_lab(0:3,1:npsi,-nmu:nmu,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    u_lab = 0.0

    ! Fluid velocities
    allocate(vfluid    (1:3,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(vfluid_old(1:3,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    vfluid     = 0.0
    vfluid_old = 0.0

    ! GR terms
    allocate(alpha      (c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(alpha_old  (c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(beta     (3,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(beta_old (3,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(phiconf    (c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(phiconf_old(c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    alpha       = 0.0
    alpha_old   = 0.0
    beta        = 0.0
    beta_old    = 0.0
    phiconf     = 0.0
    phiconf_old = 0.0

    ! Allocate main arrays
    allocate(f      (1:nflav,1:neps,1:npsi,-nmu:nmu,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(f_eq   (1:nflav,1:neps,                c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(kappa_a(1:nflav,1:neps,                c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    allocate(kappa_s(1:nflav,1:neps,                c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
    f       = 0.0
    f_eq    = 0.0
    kappa_a = 0.0
    kappa_s = 0.0

    if (heun) then
      allocate(f_old(1:nflav,1:neps,1:npsi,-nmu:nmu,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2), stat=allocstat)
      f_old = 0.0
    endif

    allocate(f_null (1:nflav,1:neps,1:npsi,-nmu:nmu                                   ), stat=allocstat)
    f_null = 0.0

    ! Jump indicator switch (store the entire array)
    allocate(jump(0:nphi,0:ntheta,0:nr), stat=allocstat)
    jump = .false.

    ! Communication send/receive buffers
    call calculate_nbuf
    allocate(sbuf(nbuf), stat=allocstat)
    sbuf = 0.0
#ifdef MPI_BT
    allocate(rbuf(nbuf), stat=allocstat)
    rbuf = 0.0
#else
    rbuf => sbuf
#endif /* MPI_BT */

  end subroutine allocate_memory

  subroutine equal_spacing(x, x_if, start, end, size)
    integer,  intent(in)  :: size
    real,     intent(out) :: x(1:size)
    real,     intent(out) :: x_if(0:size)
    real,     intent(in)  :: start
    real,     intent(in)  :: end

    integer :: i

    x_if(0) = start
    do i = 1, size
      x_if(i) = x_if(0) + (end - start) * (real(i) / real(size))
      x(i) = 0.5 * (x_if(i) + x_if(i-1))
    enddo

  end subroutine equal_spacing

  subroutine log_sphere_spacing(x, x_if, startval, endval, size)
    integer,  intent(in)  :: size
    real,     intent(out) :: x(1:size)
    real,     intent(out) :: x_if(0:size)
    real,     intent(in)  :: startval
    real,     intent(in)  :: endval

    integer :: i

    x_if(0) = 0.0

    if (size == 1) then
      x_if(1) = endval
      x(1) = (x_if(1)**4 - x_if(0)**4) / (x_if(1)**3 - x_if(0)**3) * 3.0/4.0
    else
      do i = 1, size
        x_if(i) = startval * (endval / startval)**(real(i-1) / real(size-1))
        x(i) = (x_if(i)**4 - x_if(i-1)**4) / (x_if(i)**3 - x_if(i-1)**3) * 3.0/4.0
      enddo
    endif

  end subroutine log_sphere_spacing

  subroutine calculate_rmin
    use rice_config, only: rmin, rmax, nr

    rmin = rmax**(1.0/real(nr)) - 1.0

  end subroutine calculate_rmin

  subroutine spatial_grids
    use rice_constants,       only: pi
    use rice_config,          only: rmin, rmax, thetamin, thetamax, phimin, phimax, &
                                    cartoon_grid, auto_grid, auto_rmin
    use rice_boundaries_grid, only: fill_ghosts_grid, fill_ghosts_grid_if

    real    :: thmin, thmax
    real    :: w
    integer :: b


    ! -- radial coordinates --

    if (auto_rmin) then
      print*, '(grid) automatically calculating rmin'
      call calculate_rmin
    endif

    print*, '(grid) r:      ', nr, 'points from', rmin, 'to', rmax, '(log spacing)'
    call log_sphere_spacing(r(1:nr), r_if(0:nr), rmin, rmax, nr)
    call fill_ghosts_grid(r, nr, r_if, 1)
    call fill_ghosts_grid_if(r_if, nr)


    ! -- theta coordinates --

    if ((ntheta == 1 .and. nphi == 1) .and. auto_grid) then
      call auto_grid_angles
    endif

    print*, '(grid) theta:  ', ntheta, 'points from', thetamin, 'to', thetamax
    call equal_spacing(theta(1:ntheta), theta_if(0:ntheta), thetamin, thetamax, ntheta)

    grid_2d = (ntheta > 1)

    if (grid_2d) then
      if (cartoon_grid) stop 'Cannot use cartoon grid in 2D'

      ! If grid reaches the poles, shift the cell centers so that the pole has a cell center
      if (thetamin <= epsilon(thetamin) .and. thetamax >= pi - epsilon(thetamax)) then
        print*, '(grid) shifting central theta values so that theta = 0, pi'
        thmin = theta(1)
        thmax = theta(ntheta)
        do b = 1, ntheta
          w = (theta(b) - thmin) / (thmax - thmin)
          theta(b) = theta_if(b-1) + w * (theta_if(b) - theta_if(b-1))
        enddo
      else
        stop 'Run is in 2D but thetamin and thetamax do not extend to poles'
      endif
    endif

    call fill_ghosts_grid(theta, ntheta, theta_if, 2)
    call fill_ghosts_grid_if(theta_if, ntheta)


    ! -- phi coordinates --

    print*, '(grid) phi:    ', nphi, 'points from', phimin, 'to', phimax
    call equal_spacing(phi(1:nphi), phi_if(0:nphi), phimin, phimax, nphi)
    call fill_ghosts_grid(phi, nphi, phi_if, 3)
    call fill_ghosts_grid_if(phi_if, nphi)

  end subroutine spatial_grids

  subroutine phase_grids
    use rice_config,    only: epsmin, epsmax
    use rice_constants, only: pi

    real :: psi_offset

    print*, '(grid) mu:     ', 2*nmu+1, 'at collocation points'
    call calculate_collocation_points(mu, mu_if)

    ! Offset psi so that cell center is orthogonal
    psi_offset = -0.5 * 2.0*pi / real(npsi)
    print*, '(grid) psi:    ', npsi, 'points equally spaced'
    call equal_spacing(psi, psi_if, 0.0 + psi_offset, 2.0*pi + psi_offset, npsi)

    print*, '(grid) calculating domega'
    call get_domega

    if (neps == 1) then
      eps(1) = 1.0
      eps_if(0) = 0.0
      eps_if(1) = 3.0**(1.0/3.0)
    else
      print*, '(grid) epsilon:', neps, 'points from', epsmin, 'to', epsmax
      call log_sphere_spacing(eps, eps_if, epsmin, epsmax, neps)
    endif

    print*, '(grid) calculating phase space volume'
    call get_epsvol

  end subroutine phase_grids

  subroutine aux_grids
    use rice_config, only: cartoon_grid, check_grid

    print*, '(grid) calculating interface areas'
    call get_area

    print*, '(grid) calculating cell volumes'
    call get_volume

    if (cartoon_grid) then
      print*, '(grid) overriding interface areas and volumes with equatorial values'
      call cartoon_override
    else if (.not. external_coordinates) then
      print*, '(grid) calculating mid-cell theta angle'
      ! Overwrites theta_if, must do after area
      call get_dtheta_mid
    endif

    if (.not. external_coordinates) then
      print*, '(grid) checking interface areas'
      call check_area
    else
      print*, '(grid) using external coordinates, not checking cell area'
    endif

    print*, '(grid) calculating cell ratio'
    ! Must calculate after angles are set
    call get_cell_ratio

    if (check_grid) then
      print*, '(grid) checking cell aspect ratio'
      call check_cell_aspect_ratio
    endif

    print*, '(grid) calculating comoving velocities for all bins in phase space'
    call get_ucom

    print*, '(grid) precomputing trig values for psi'
    call precompute_trig

    print*, '(grid) precomputing differentials'
    call precompute_differentials

    if (cartoon_grid) then
      print*, '(grid) precomputing cartoon rotation weights'
      call precompute_cartoon_weights
    endif

    print*, '(grid) setting f_null to zero'
    f_null(:,:,:,:) = 0.0

  end subroutine aux_grids

  subroutine calculate_collocation_points(mu, mu_if)
    use rice_constants,    only: pi
    use rice_gausslobatto, only: collocation_points
    use rice_config,       only: nmu, uniform_mu

    real,     intent(out) :: mu(-nmu:nmu)
    real,     intent(out) :: mu_if(-nmu-1:nmu)

    real :: theta_mom(-nmu:nmu)
    real :: theta_mom_if(-nmu-1:nmu)

    real :: wtmu(-nmu:nmu)

    integer :: i

    ! Give more resolution at mu=0, but worse accuracy overall
    if (uniform_mu) then
      call equal_spacing(theta_mom, theta_mom_if, 0.0, pi, 2*nmu+1)

      do i = -nmu, nmu
        mu(i)    = - cos(theta_mom(i))
      enddo

      do i = -nmu-1, nmu
        mu_if(i) = - cos(theta_mom_if(i))
      enddo

      ! Point ends to -1 and 1
      mu(-nmu) = -1.0
      mu( nmu) =  1.0

    else
      call collocation_points(2*nmu+1, mu(-nmu:nmu), wtmu(-nmu:nmu))

      mu_if(-nmu-1) = -1.0
      do i = -nmu, nmu
        mu_if(i) = mu_if(i-1) + wtmu(i)
      enddo
    endif

    print*, '(grid) adjusting mu to be symmetric'
    mu(0) = 0.0
    do i = 1, nmu
      mu(i) = -mu(-i)
    enddo
    do i = 0, nmu
      mu_if(i) = -mu_if(-(i+1))
    enddo

  end subroutine calculate_collocation_points

  subroutine get_area
    use rice_config, only: closed_boundaries, square_inner
    use rice_utils,  only: cross_product, magnitude, cartesian_coordinate

    integer :: a, b, c

    real :: dr
    real :: dtheta
    real :: dphi
    real :: l_phi_tl
    real :: l_phi_tr
    real :: l_phi_rl
    real :: l_phi_rr
    real :: l_theta_rl
    real :: l_theta_rr
    real :: l_theta

    real :: v1(3)
    real :: v2(3)
    real :: v3(3)
    real :: v4(3)

    real :: area_numeric

    ! Angles are measured between vertices
    ! Lengths are measured between vertices

    do a = 0, nr
      do b = 1, ntheta
        do c = 1, nphi
          ! Exact area
          dphi = phi_if(c) - phi_if(c-1)
          dtheta = theta_if(b) - theta_if(b-1)
          l_phi_tl = 2.0 * r_if(a) * sin(0.5*dphi) * sin(theta_if(b-1))
          l_phi_tr = 2.0 * r_if(a) * sin(0.5*dphi) * sin(theta_if(b  ))
          l_theta = 2.0 * r_if(a) * sin(0.5*dtheta)
          area_r(c,b,a) = 0.5 * (l_phi_tl + l_phi_tr) * l_theta * cos(0.5*dtheta)

          ! Numeric area using triangles
          v1 = cartesian_coordinate(r_if(a), theta_if(b-1), phi_if(c-1))
          v2 = cartesian_coordinate(r_if(a), theta_if(b),   phi_if(c-1))
          v3 = cartesian_coordinate(r_if(a), theta_if(b-1), phi_if(c)  )
          v4 = cartesian_coordinate(r_if(a), theta_if(b  ), phi_if(c)  )
          area_numeric = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

          area_r(c,b,a) = area_numeric

        enddo
      enddo
    enddo

    ! Make area at center exactly zero
    area_r(:,:,0) = 0.0

    ! Disable outflow for conservation testing
    if (closed_boundaries) then
      area_r(:,:,nr) = 0.0
    endif

    do a = 1, nr
      do b = 0, ntheta
        do c = 1, nphi
          ! Exact area
          dr = r_if(a) - r_if(a-1)
          dphi = phi_if(c) - phi_if(c-1)
          l_phi_rl = 2.0 * r_if(a-1) * sin(0.5*dphi) * sin(theta_if(b))
          l_phi_rr = 2.0 * r_if(a  ) * sin(0.5*dphi) * sin(theta_if(b))
          area_theta(c,b,a) = 0.5 * (l_phi_rl + l_phi_rr) * dr * cos(0.5*dphi) ! incorrect

          ! Numeric area using triangles
          v1 = cartesian_coordinate(r_if(a-1), theta_if(b), phi_if(c-1))
          v2 = cartesian_coordinate(r_if(a-1), theta_if(b), phi_if(c)  )
          v3 = cartesian_coordinate(r_if(a)  , theta_if(b), phi_if(c-1))
          v4 = cartesian_coordinate(r_if(a)  , theta_if(b), phi_if(c)  )

          area_numeric = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

          area_theta(c,b,a) = area_numeric

        enddo
      enddo
    enddo

    ! Make areas at poles exactly zero
    if (grid_2d) then
      area_theta(:,0,     :) = 0.0
      area_theta(:,ntheta,:) = 0.0
    endif

    do a = 1, nr
      do b = 1, ntheta
        do c = 0, nphi
          ! Exact area
          dr = r_if(a) - r_if(a-1)
          dtheta = theta_if(b) - theta_if(b-1)
          l_theta_rl = 2.0 * r_if(a-1) * sin(0.5*dtheta)
          l_theta_rr = 2.0 * r_if(a  ) * sin(0.5*dtheta)
          area_phi(c,b,a) = 0.5 * (l_theta_rl + l_theta_rr) * dr * cos(0.5*dtheta)

          ! Numeric area using triangles
          v1 = cartesian_coordinate(r_if(a-1), theta_if(b-1), phi_if(c))
          v2 = cartesian_coordinate(r_if(a-1), theta_if(b),   phi_if(c))
          v3 = cartesian_coordinate(r_if(a)  , theta_if(b-1), phi_if(c))
          v4 = cartesian_coordinate(r_if(a)  , theta_if(b  ), phi_if(c))
          area_numeric = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

          area_phi(c,b,a) = area_numeric

        enddo
      enddo
    enddo

    ! Hack innermost zone to be square
    if (square_inner) then
      area_r(:,:,0) = area_r(:,:,1)
    endif

  end subroutine get_area

  subroutine check_area
    use rice_config, only: nr, ntheta, nphi

    real :: diff
    real :: area_lat

    integer :: a, b, c

    if (ntheta == 1) then

      do a = 2, nr-1
        do b = 1, ntheta
          do c = 1, nphi
            diff = area_r(c,b,a) - area_r(c,b,a-1)
            area_lat = 0.0
            area_lat = area_lat + sin(theta(b)-theta_if(b-1)) * area_theta(c,b-1,a)
            area_lat = area_lat + sin(theta_if(b)-theta(b)  ) * area_theta(c,b,a)
            area_lat = area_lat + sin(phi(c)-phi_if(c-1))     * area_phi(c-1,b,a)
            area_lat = area_lat + sin(phi_if(c)-phi(c)  )     * area_phi(c,b,a)

            if (abs(1.0 - diff/area_lat) > 1.e-13) then
              print*, '(grid) cell area check failed'
              stop
            endif

          enddo
        enddo
      enddo

    endif

  end subroutine check_area

  subroutine get_volume
    use rice_utils,           only: determinant, cartesian_coordinate
    use rice_boundaries_grid, only: fill_ghosts_scalar_1
    integer :: a, b, c

    real :: dr
    real :: dtheta
    real :: dphi

    real :: v0(3)
    real :: v1(3)
    real :: v2(3)
    real :: v3(3)
    real :: v4(3)
    real :: v5(3)
    real :: v6(3)
    real :: v7(3)

    real :: exact

    real :: matrix(3,3)

    real :: height

    exact = 0.0

    do a = 1, nr
      do b = 1, ntheta
        do c = 1, nphi

          ! Analytic volume calculation - there is likely an error here
          dr = r_if(a) - r_if(a-1)
          dtheta = theta_if(b) - theta_if(b-1)
          dphi = phi_if(c) - phi_if(c-1)
          ! height = dr * sqrt(0.5*(cos(dtheta) + cos(dphi)))
          height = dr * cos(0.5*dtheta) * cos(0.5*dphi)
          exact = exact + 0.5 * (area_r(c,b,a-1) + area_r(c,b,a)) * height

          ! The most bulletproof way is to divide the general hexahedron into tetrahedra
          v0 = cartesian_coordinate(r_if(a-1), theta_if(b-1), phi_if(c-1))
          v1 = cartesian_coordinate(r_if(a-1), theta_if(b  ), phi_if(c-1))
          v2 = cartesian_coordinate(r_if(a-1), theta_if(b-1), phi_if(c  ))
          v3 = cartesian_coordinate(r_if(a-1), theta_if(b  ), phi_if(c  ))
          v4 = cartesian_coordinate(r_if(a  ), theta_if(b-1), phi_if(c-1))
          v5 = cartesian_coordinate(r_if(a  ), theta_if(b  ), phi_if(c-1))
          v6 = cartesian_coordinate(r_if(a  ), theta_if(b-1), phi_if(c  ))
          v7 = cartesian_coordinate(r_if(a  ), theta_if(b  ), phi_if(c  ))

          volume(c,b,a) = 0.0

          matrix(:,1) = v7 - v0
          matrix(:,2) = v1 - v0
          matrix(:,3) = v3 - v5
          volume(c,b,a) = volume(c,b,a) + determinant(matrix)

          matrix(:,1) = v7 - v0
          matrix(:,2) = v4 - v0
          matrix(:,3) = v5 - v6
          volume(c,b,a) = volume(c,b,a) + determinant(matrix)

          matrix(:,1) = v7 - v0
          matrix(:,2) = v2 - v0
          matrix(:,3) = v6 - v3
          volume(c,b,a) = volume(c,b,a) + determinant(matrix)

          volume(c,b,a) = volume(c,b,a) / 6.0

        enddo
      enddo
    enddo

    print*, '(grid) exact volume / brute force volume = ', exact / sum(volume)

    ! Get ghost zone volumes
    do a = 1, nr
      do b = 1, ntheta
        call fill_ghosts_scalar_1(volume(:,b,a), nphi, 3)
      enddo
    enddo
    do a = 1, nr
      do c = 0, nphi+1 ! fill corners
        call fill_ghosts_scalar_1(volume(c,:,a), ntheta, 2)
      enddo
    enddo
    do b = 0, ntheta+1 ! fill corners
      do c = 0, nphi+1 ! fill corners
        call fill_ghosts_scalar_1(volume(c,b,:), nr, 1)
      enddo
    enddo

  end subroutine get_volume

  subroutine get_domega
    integer :: j, k
    omega = 0.0
    do j = -nmu, nmu
      do k = 1, npsi
        domega(k,j) = (psi_if(k) - psi_if(k-1)) * (mu_if(j) - mu_if(j-1))
        omega = omega + domega(k,j)
      enddo
    enddo
    domega_scaled = domega / omega
    domega_inv = 1.0 / domega
  end subroutine get_domega

  subroutine get_epsvol
    integer :: i

    do i = 1, neps
      epsvol(i) = (1.0/3.0) * (eps_if(i)**3 - eps_if(i-1)**3)
    enddo

    epsvol_inv = 1.0 / epsvol

  end subroutine get_epsvol

  subroutine get_ucom
    integer :: j, k
    real    :: sintheta

    do j = -nmu, nmu
      sintheta = sqrt(1.0 - mu(j)**2)
      do k = 1, npsi
        u_com(0,k,j) = 1.0
        u_com(1,k,j) = mu(j)
        u_com(2,k,j) = sintheta * cos(psi(k))
        u_com(3,k,j) = sintheta * sin(psi(k))

        if (abs(u_com(1,k,j)) <= epsilon(u_com(1,k,j))) u_com(1,k,j) = 0.0
        if (abs(u_com(2,k,j)) <= epsilon(u_com(2,k,j))) u_com(2,k,j) = 0.0
        if (abs(u_com(3,k,j)) <= epsilon(u_com(3,k,j))) u_com(3,k,j) = 0.0

      enddo
    enddo

  end subroutine get_ucom

  subroutine precompute_trig
    integer :: i
    do i = 1, npsi
       cpsi(i) = cos(psi (i))
       spsi(i) = sin(psi (i))
    end do

  end subroutine precompute_trig

  subroutine get_dtheta_mid
    use rice_utils, only: magnitude, cartesian_coordinate
    use rice_boundaries_grid, only: fill_ghosts_grid, fill_ghosts_grid_if

    integer, parameter :: qp = selected_real_kind(33, 4931)
    integer, parameter :: dp = selected_real_kind(15, 307)

    real :: v1(3), v2(3), v3(3), v4(3)
    real :: vmid_top(3), vmid_bottom(3), vmid(3)

    real, allocatable :: theta_if_new(:)

    real :: z(3)

    integer :: b

    allocate(theta_if_new(0:ntheta))

    z(1) = 0.0
    z(2) = 0.0
    z(3) = 1.0

    do b = 1, ntheta
      ! Assume evenly spaced phi grid
      v1 = cartesian_coordinate(1.0, theta_if(b-1), phi_if(0))
      v2 = cartesian_coordinate(1.0, theta_if(b),   phi_if(0))
      v3 = cartesian_coordinate(1.0, theta_if(b-1), phi_if(1))
      v4 = cartesian_coordinate(1.0, theta_if(b  ), phi_if(1))

      vmid_top = 0.5 * (v2 + v4)
      vmid_bottom = 0.5 * (v1 + v3)
      vmid = 0.5 * (vmid_top + vmid_bottom)

      ! Calculate acos() using quad precision for consistency across different platforms
      theta_if_new(b-1) = real(acos( real(dot_product(z, vmid_bottom) / (magnitude(z)*magnitude(vmid_bottom)), kind=qp) ), kind=dp)
      theta_if_new(b)   = real(acos( real(dot_product(z, vmid_top)    / (magnitude(z)*magnitude(vmid_top)),    kind=qp) ), kind=dp)
    enddo

    ! Write after all done, otherwise will interfere with loop
    theta_if(0:ntheta) = theta_if_new(0:ntheta)

    deallocate(theta_if_new)

    call fill_ghosts_grid(theta, ntheta, theta_if, 2, cover_poles=grid_2d)
    call fill_ghosts_grid_if(theta_if, ntheta)

  end subroutine get_dtheta_mid

  subroutine get_cell_ratio

    real :: radial_area
    integer :: c,b,a

    do a = 1, nr
      do b = 1, ntheta
        do c = 1, nphi
          radial_area = area_r(c,b,a-1)
          radial_area = radial_area + area_theta(c,b-1,a) * sin(theta(b)    - theta_if(b-1))
          radial_area = radial_area + area_theta(c,b,  a) * sin(theta_if(b) - theta(b)     )
          radial_area = radial_area + area_phi(c-1,b,a)   * sin(phi(c)      - phi_if(c-1)  )
          radial_area = radial_area + area_phi(c,  b,a)   * sin(phi_if(c-1) - phi_if(c)    )

          ! How much of the radial flux is through angular interfaces
          cell_ratio(c,b,a) = 1.0 - max(0.0, area_r(c,b,a-1) / radial_area)
        enddo
      enddo
    enddo

  end subroutine get_cell_ratio

  subroutine precompute_differentials
    integer :: i, j, k

    do j = -nmu, nmu
      do k = 1, npsi
        do i = 1, neps
          eps2domegadeps   (i,k,j) = domega(k,j) * epsvol(i)
          eps3domegadeps   (i,k,j) = eps(i) * eps2domegadeps(i,k,j)
          mueps3domegadeps (i,k,j) = mu(j) * eps3domegadeps(i,k,j)
          mu2eps3domegadeps(i,k,j) = mu(j) * mueps3domegadeps(i,k,j)
        enddo
      enddo
    enddo
  end subroutine precompute_differentials

  subroutine cartoon_override
    use rice_config, only: ntheta, thetamin, thetamax, closed_boundaries
    use rice_utils,  only: cross_product, magnitude, determinant, cartesian_coordinate

    real :: theta_l
    real :: theta_r

    real :: v0(3)
    real :: v1(3)
    real :: v2(3)
    real :: v3(3)
    real :: v4(3)
    real :: v5(3)
    real :: v6(3)
    real :: v7(3)

    real :: matrix(3,3)

    integer :: a, b, c

    integer :: allocstat

    ! Save original volumes and areas
    allocate(volume_spherical    (1:nphi,1:ntheta,1:nr), stat=allocstat)
    allocate(area_r_spherical    (1:nphi,1:ntheta,0:nr), stat=allocstat)
    allocate(area_theta_spherical(1:nphi,0:ntheta,1:nr), stat=allocstat)
    allocate(area_phi_spherical  (0:nphi,1:ntheta,1:nr), stat=allocstat)

    volume_spherical     = volume(1:nphi,1:ntheta,1:nr)
    area_r_spherical     = area_r
    area_theta_spherical = area_theta
    area_phi_spherical   = area_phi

    if (ntheta == 1) then
      print*, 'Error: Cartoon grid needs to be 2D'
      stop
    endif

    if (nphi > 1) then
      print*, 'Error: Cartoon grid cannot run in 3D'
      stop
    endif

    theta_l = 0.5*(thetamin + thetamax) - 0.5 * (thetamax - thetamin) / real(ntheta)
    theta_r = 0.5*(thetamin + thetamax) + 0.5 * (thetamax - thetamin) / real(ntheta)

    c = 1

    do a = 0, nr
      do b = 1, ntheta
        ! Numeric area using triangles
        v1 = cartesian_coordinate(r_if(a), theta_l, phi_if(c-1))
        v2 = cartesian_coordinate(r_if(a), theta_r, phi_if(c-1))
        v3 = cartesian_coordinate(r_if(a), theta_l, phi_if(c)  )
        v4 = cartesian_coordinate(r_if(a), theta_r, phi_if(c)  )

        area_r(c,b,a) = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

      enddo
    enddo

    ! Make area at center exactly zero
    area_r(:,:,0) = 0.0

    ! Disable outflow for conservation testing
    if (closed_boundaries) then
      area_r(:,:,nr) = 0.0
    endif

    do a = 1, nr
      do b = 0, ntheta
        ! Numeric area using triangles
        v1 = cartesian_coordinate(r_if(a-1), theta_r, phi_if(c-1))
        v2 = cartesian_coordinate(r_if(a-1), theta_r, phi_if(c)  )
        v3 = cartesian_coordinate(r_if(a)  , theta_r, phi_if(c-1))
        v4 = cartesian_coordinate(r_if(a)  , theta_r, phi_if(c)  )

        area_theta(c,b,a) = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

      enddo
    enddo

    ! Make pole area zero
    ! area_theta(:,0     ,:) = 0.0
    ! area_theta(:,ntheta,:) = 0.0

    do a = 1, nr
      do b = 1, ntheta
        ! Numeric area using triangles
        v1 = cartesian_coordinate(r_if(a-1), theta_l, phi_if(c))
        v2 = cartesian_coordinate(r_if(a-1), theta_r, phi_if(c))
        v3 = cartesian_coordinate(r_if(a)  , theta_l, phi_if(c))
        v4 = cartesian_coordinate(r_if(a)  , theta_r, phi_if(c))

        area_phi(c,b,a) = 0.5 * (magnitude(cross_product(v2-v1,v3-v1)) + magnitude(cross_product(v2-v4,v3-v4)))

      enddo
    enddo

    do a = 1, nr
      do b = 1, ntheta

        v0 = cartesian_coordinate(r_if(a-1), theta_l, phi_if(c-1))
        v1 = cartesian_coordinate(r_if(a-1), theta_r, phi_if(c-1))
        v2 = cartesian_coordinate(r_if(a-1), theta_l, phi_if(c  ))
        v3 = cartesian_coordinate(r_if(a-1), theta_r, phi_if(c  ))
        v4 = cartesian_coordinate(r_if(a  ), theta_l, phi_if(c-1))
        v5 = cartesian_coordinate(r_if(a  ), theta_r, phi_if(c-1))
        v6 = cartesian_coordinate(r_if(a  ), theta_l, phi_if(c  ))
        v7 = cartesian_coordinate(r_if(a  ), theta_r, phi_if(c  ))

        volume(c,b,a) = 0.0

        matrix(:,1) = v7 - v0
        matrix(:,2) = v1 - v0
        matrix(:,3) = v3 - v5
        volume(c,b,a) = volume(c,b,a) + determinant(matrix)

        matrix(:,1) = v7 - v0
        matrix(:,2) = v4 - v0
        matrix(:,3) = v5 - v6
        volume(c,b,a) = volume(c,b,a) + determinant(matrix)

        matrix(:,1) = v7 - v0
        matrix(:,2) = v2 - v0
        matrix(:,3) = v6 - v3
        volume(c,b,a) = volume(c,b,a) + determinant(matrix)

        volume(c,b,a) = volume(c,b,a) / 6.0

      enddo
    enddo

  end subroutine cartoon_override

  subroutine precompute_cartoon_weights
    use rice_config,          only: ntheta
    use rice_transform_utils, only: cartoon_rotation_angle, rotate_vector_cartoon, &
                                    get_nearest_index, get_nearest_index_reverse

    real :: beta, alpha_angle(4)

    real :: u_new(0:3)
    real :: idenom, sinpsii, cospsii
    real :: sinpsia, sinpsib, cospsia, cospsib

    real :: wpsia, wpsib, itotal

    integer :: allocstat
    integer :: b, ia, k
    integer :: ka, kb

    ! For each theta
    ! For 4 alpha angles
    ! wpsia for each psi
    allocate(cartoon_weight (npsi,4,ntheta), stat=allocstat)
    allocate(cartoon_index(2,npsi,4,ntheta), stat=allocstat)

    alpha_angle(1) = -2.0 * (phi_if(1) - phi_if(0))
    alpha_angle(2) = -1.0 * (phi_if(1) - phi_if(0))
    alpha_angle(3) =  1.0 * (phi_if(1) - phi_if(0))
    alpha_angle(4) =  2.0 * (phi_if(1) - phi_if(0))

    do b = 1, ntheta
      do ia = 1, 4
        beta = cartoon_rotation_angle(theta(b), alpha_angle(ia))

        do k = 1, npsi
          call rotate_vector_cartoon(u_com(:,k,0), u_new, sin(beta), cos(beta))
          idenom = 1.0 / sqrt(u_new(2)**2 + u_new(3)**2)
          cospsii = u_new(2) * idenom
          sinpsii = u_new(3) * idenom

          if (sinpsii > 0.5*sqrt(2.0)) then
            call get_nearest_index_reverse(cospsii, cpsi, npsi, k, 0.0, ka, kb)
          elseif (sinpsii < -0.5*sqrt(2.0)) then
            call get_nearest_index(cospsii, cpsi, npsi, k, 0.0, ka, kb)
          elseif (cospsii > 0.5*sqrt(2.0)) then
            call get_nearest_index(sinpsii, spsi, npsi, k, 0.0, ka, kb)
          else ! (cospsii < -0.5*sqrt(2.0))
            call get_nearest_index_reverse(sinpsii, spsi, npsi, k, 0.0, ka, kb)
          endif

          ! Only psi changes for this rotation
          cospsia = cpsi(ka)
          cospsib = cpsi(kb)
          sinpsia = spsi(ka)
          sinpsib = spsi(kb)

          idenom = 1.0 / (cospsib*sinpsia-cospsia*sinpsib)

          wpsib = (cospsii*sinpsia - cospsia*sinpsii) * idenom
          wpsia = (sinpsii*cospsib - sinpsib*cospsii) * idenom

          itotal = 1.0 / (wpsia + wpsib)
          wpsib = wpsib * itotal
          wpsia = wpsia * itotal

          cartoon_weight(k,ia,b) = wpsia
          cartoon_index(1,k,ia,b) = ka
          cartoon_index(2,k,ia,b) = kb
        enddo
      enddo
    enddo


  end subroutine precompute_cartoon_weights

  subroutine check_cell_aspect_ratio
    use rice_constants, only: pi
    use rice_config, only: nr, nmu, npsi

    real :: r_theta, r_phi
    integer :: a, b, c

    logical :: error

    real, parameter :: squareness = 2.0

    a = nr-1 ! using nr will cause problems when using a closed grid
    b = max(1, ntheta / 2)
    c = 1

    r_theta = area_r(c,b,a) / area_theta(c,b,a)
    r_phi   = area_r(c,b,a) / area_phi(c,b,a)

    error = .false.

    if (r_theta > squareness) then
      print*, 'Aspect ratio error: Cells too wide, increase ntheta by factor', r_theta
      error = .true.
    endif

    if (r_phi > squareness) then
      print*, 'Aspect ratio error: Cells too wide, increase nphi by factor', r_phi
      error = .true.
    endif

    if (r_theta < 1.0 / squareness) then
      print*, 'Aspect ratio error: Cells too long, decrease ntheta by factor', 1.0 / r_theta
      error = .true.
    endif

    if (r_phi < 1.0 / squareness) then
      print*, 'Aspect ratio error: Cells too long, decrease nphi by factor', 1.0 / r_phi
      error = .true.
    endif

    if (r_theta > squareness * r_phi) then
      print*, 'Aspect ratio error: Cells have bad aspect ratio, increase nphi by factor', r_theta / r_phi
      error = .true.
    elseif (r_phi > squareness * r_theta) then
      print*, 'Aspect ratio error: Cells have bad aspect ratio, increase ntheta by factor', r_phi / r_theta
      error = .true.
    endif

    ! Check that dtheta and dphi are matched with nmu
    r_theta = (theta_if(1) - theta_if(0)) / (pi / real(2*nmu+1))
    r_phi = (phi_if(1) - phi_if(0)) / (pi / real(2*nmu+1))

    if (max(r_theta, 1.0/r_theta) > 2.0) then
      print*, '(grid) dtheta is highly mismatched with nmu, solution may be inaccurate. Stopping.'
      error = .true.
    endif

    if (max(r_phi, 1.0/r_phi) > 2.0) then
      print*, '(grid) dphi is highly mismatched with nmu, solution may be inaccurate. Stopping.'
      error = .true.
    endif

    ! Check that dphi is matched with npsi
    r_phi = (phi_if(1) - phi_if(0)) / (2.0*pi / real(npsi))
    if (max(r_phi, 1.0/r_phi) > 2.0) then
      print*, '(grid) dphi is highly mismatched with npsi, solution may be inaccurate. Stopping.'
      error = .true.
    endif

    if (error) stop

  end subroutine check_cell_aspect_ratio

  subroutine auto_grid_angles
    use rice_constants, only: pi
    use rice_config,    only: thetamin, thetamax, phimin, phimax

    real    :: dangle

    print*, '(grid) Automatically setting dtheta and dphi (ignoring config values)'

    dangle = (r_if(nr) - r_if(nr-1)) / r_if(nr-1)

    ! Adjust so that limiting direction is still 1 for efficiency
    dangle = 1.5 * dangle

    thetamin = 0.5*pi - 0.5*dangle
    thetamax = 0.5*pi + 0.5*dangle
    phimin = 0.0
    phimax = dangle

    print*, '(grid)   thetamin =', thetamin
    print*, '(grid)   thetamax =', thetamax
    print*, '(grid)   phimin   =', phimin
    print*, '(grid)   phimin   =', phimax

  end subroutine auto_grid_angles

  subroutine calculate_nbuf
    use rice_config, only: nmu, npsi, neps, nflav

    integer :: n_space ! memnory required for one point in space

    n_space = neps * nflav                    & ! f_eq
            + neps * nflav                    & ! kappa_a
            + neps * nflav                    & ! kappa_s
            + 3                               & ! vfluid
            + 1                               & ! alpha
            + 3                               & ! beta
            + 1                               & ! phiconf
            + (2*nmu+1) * npsi * neps * nflav   ! f

    nbuf = n_space * max(max((a_e-a_s+5)*(b_e-b_s+5), (a_e-a_s+5)*(c_e-c_s+5)), (b_e-b_s+5)*(c_e-c_s+5))
    nbuf = nbuf * 2 ! 2 zones deep

  end subroutine calculate_nbuf

end module rice_grid
