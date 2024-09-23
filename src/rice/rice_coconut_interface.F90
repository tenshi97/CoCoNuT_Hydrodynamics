module rice_coconut_interface
  use phycon, only: pc_pi, pc_geog, pc_kmev, pc_mb, pc_cl, pc_meverg, pc_gc, &
                    wc_s0, wc_ga2, wc_mb, wc_mn, wc_me, wc_mp, wc_mq, &
                    wc_cv, wc_hc, wc_ca
  use rice_config, only: use_boltzmann, nr

  implicit none

  ! Constant used for calculating f_eq
  ! real, parameter :: bnorm = 4.0*pc_pi*pc_cl / (2.0 * pc_pi * wc_hc)**3
  real, parameter :: bnorm = 1.0

  real, parameter :: wc_mnc=0.5*(wc_mn+wc_mp)

  real, parameter :: third = 1.0 / 3.0

  ! Density threshold for enabling bremsstrahlung term
  ! 1.e9 is suitable, but for performance profiling,
  ! set this to 0.0 so that the load is indicative of post-bounce
  ! conditions.
  real, parameter :: denbrems = 0.0

  real :: dt_rad

contains

  subroutine boltzmann_mpi_setup(nprocs_in, ntheta_in, nphi_in)
    use rice_config, only: mpi_nr, mpi_ntheta, mpi_nphi
    use rice_mpi,    only: setup_mpi
    use mpi_domains, only: optimal_decomposition

    integer, intent(in) :: nprocs_in
    integer, intent(in) :: ntheta_in
    integer, intent(in) :: nphi_in

    ! Call the same domain decomposition routine as hydro
    if (nphi_in > 1) then
      call optimal_decomposition(ntheta_in, nphi_in, nprocs_in, mpi_ntheta, mpi_nphi)
    else
      mpi_ntheta = nprocs_in
      mpi_nphi = 1
    endif
    mpi_nr = 1

    ! Setup MPI
    call setup_mpi

  end subroutine boltzmann_mpi_setup

  subroutine allocate_boltzmann_solver(nr_in, ntheta_in, nphi_in, nflav_in, qy_s, qy_e, qz_s, qz_e)
    use rice_config, only: nr, ntheta, nphi, neps, nflav
    use rice_grid,   only: allocate_memory, set_domain, &
                           a_s, a_e, na, &
                           b_s, b_e, nb, &
                           c_s, c_e, nc
    use rice_mpi,    only: calculate_ncyc

    integer,  intent(in)  :: nr_in
    integer,  intent(in)  :: ntheta_in
    integer,  intent(in)  :: nphi_in
    integer,  intent(in)  :: nflav_in
    integer,  intent(in)  :: qy_s, qy_e, qz_s, qz_e

    ! Array sizes
    nr = nr_in
    ntheta = ntheta_in
    nphi = nphi_in

    ! Copy domain from hydro
    a_s = 1
    a_e = nr
    na = nr

    b_s = qy_s
    b_e = qy_e
    nb = b_e - b_s + 1

    c_s = qz_s
    c_e = qz_e
    nc = c_e - c_s + 1

    ! Allocate grids
    call allocate_memory

    ! Calculate number of cycles based on ray width
    call calculate_ncyc

  end subroutine allocate_boltzmann_solver

  subroutine init_boltzmann_solver(r_in, theta_in, phi_in, &
                                   r_if_in, theta_if_in, phi_if_in, clight_in)
    use rice_constants,       only: pi
    use rice_config,          only: nr, ntheta, nphi, neps, nmu, npsi, clight, &
                                    epsmin, epsmax, nmu, thetamin, thetamax, phimin, phimax, &
                                    check_grid
    use rice_grid,            only: phase_grids, aux_grids, &
                                    r, theta, phi, r_if, theta_if, phi_if, &
                                    f, time, mu, vfluid_old, alpha_old, beta_old, phiconf_old, &
                                    external_coordinates, grid_2d
    use rice_output,          only: read_output
    use rice_dt,              only: get_dt
    use rice_tune,            only: tune_areas
    use rice_boundaries_grid, only: fill_ghosts_grid, fill_ghosts_grid_if
    real,     intent(in)  :: r_in(-1:nr+2)
    real,     intent(in)  :: theta_in(-1:ntheta+2)
    real,     intent(in)  :: phi_in(-1:nphi+2)
    real,     intent(in)  :: r_if_in(-1:nr+1)
    real,     intent(in)  :: theta_if_in(-1:ntheta+1)
    real,     intent(in)  :: phi_if_in(-1:nphi+1)
    real,     intent(in)  :: clight_in

    real :: wedge_angle

    print*, '(bt) Initializing Boltzmann solver'

    ! Config options
    epsmin = 4.0
    epsmax = 240.0

    ! Grid coordinates
    ! Inherit the ghost coordinates from hydro
    r     (-1:nr+2)     = r_in     (-1:nr+2)
    theta (-1:ntheta+2) = theta_in (-1:ntheta+2)
    phi   (-1:nphi+2)   = phi_in   (-1:nphi+2)

    r_if     (-1:nr+1)     = r_if_in     (-1:nr+1)
    theta_if (-1:ntheta+1) = theta_if_in (-1:ntheta+1)
    phi_if   (-1:nphi+1)   = phi_if_in   (-1:nphi+1)

    ! Tell the BT solver to use the specified coordinates
    ! Prevents interfaces from being adjusted,
    ! unless the dimension only has one cell
    external_coordinates = .true.

    grid_2d = (ntheta > 1)

    call phase_grids

    ! Set theta and phi wedge sizes
    wedge_angle = acos(mu(-nmu)) - acos(mu(-nmu+1))

    if (ntheta == 1) then
      thetamax = 0.5*pi + 0.5*wedge_angle
      thetamin = 0.5*pi - 0.5*wedge_angle
      theta(1) = 0.5*(thetamin + thetamax)
      theta_if(0) = thetamin
      theta_if(1) = thetamax
      call fill_ghosts_grid(theta, ntheta, theta_if, 2)
      call fill_ghosts_grid_if(theta_if, ntheta)
    endif

    if (nphi == 1) then
      phimax = wedge_angle
      phimin = 0.0
      phi = 0.5*(phimin + phimax)
      phi_if(0) = phimin
      phi_if(1) = phimax
      call fill_ghosts_grid(phi, nphi, phi_if, 3)
      call fill_ghosts_grid_if(phi_if, nphi)
    endif

    check_grid = .false.

    call aux_grids

    ! call tune_areas(0.962272232181095)

    ! Set speed of light
    clight = clight_in

    ! Get radiation dt
    call get_dt(dt_rad)

    if (time == 0.0) then ! If this is a new start
      ! Initialise f to zero
      print*, '(bt) Initializing f to zero'
      f = 0.0

      print*, '(bt) Initializing vfluid_old to zero'
      vfluid_old = 0.0

      print*, '(bt) Initializing alpha_old to one'
      alpha_old = 1.0

      print*, '(bt) Initializing beta_old to zero'
      beta_old = 0.0

      print*, '(bt) Initializing phiconf_old to one'
      phiconf_old = 1.0
    endif

  end subroutine init_boltzmann_solver

  subroutine boltzmann_step(t_start, dt_hydro, rho, t, xnnu, nx, &
                            v_1, v_2, v_3, &
                            phi, alpha_in, beta_1, beta_2, beta_3, &
                            qye, qen, qmo_r, &
                            dnu, enu, fnu, pnu, enu_lab, fnu_lab)
    use rice_config,    only: neps, nmu, npsi, nflav, &
                              metric_contribution, enable_gr
    use rice_step,      only: step, time_boost
    use rice_output,    only: write_output
    use rice_grid,      only: f, vfluid, vfluid_old, alpha, alpha_old, &
                              beta, beta_old, phiconf, phiconf_old, time, &
                              a_s, a_e, b_s, b_e, c_s, c_e

    ! CoCoNuT indexing is the other way around
    real,    intent(in)    :: t_start
    real,    intent(in)    :: dt_hydro
    integer, intent(in)    :: nx
    real,    intent(in)    :: rho     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: t       (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: xnnu    (a_s:a_e,b_s:b_e,c_s:c_e,1:nx)
    real,    intent(in)    :: v_1     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: v_2     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: v_3     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: phi     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: alpha_in(a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: beta_1  (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: beta_2  (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(in)    :: beta_3  (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(inout) :: qye     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(inout) :: qen     (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(inout) :: qmo_r   (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(out)   :: dnu     (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real,    intent(out)   :: enu     (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real,    intent(out)   :: fnu     (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real,    intent(out)   :: pnu     (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real,    intent(out)   :: enu_lab (a_s:a_e,b_s:b_e,c_s:c_e)
    real,    intent(out)   :: fnu_lab (a_s:a_e,b_s:b_e,c_s:c_e)

    real                   :: qmo   (3,a_s:a_e,b_s:b_e,c_s:c_e)

    real :: factor
    real :: insteps
    real :: dt

    real :: const
    real :: constc

    integer :: a, b, c
    integer :: i, j, k, l

    integer :: nsteps
    integer :: istep

    ! Update time step, put it back into exact sync
    time = t_start + dt_hydro

    ! Number of substeps
    nsteps = ceiling(dt_hydro/dt_rad)

    ! Decrease dt so that it is an integer multiple of dt_hydro
    dt = dt_hydro / real(nsteps)

    ! Convert velocity to CGS units
    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e
          ! Order is switched
          factor = pc_cl / phi(a,b,c)**2
          vfluid(1,c,b,a) = v_1(a,b,c) * factor
          vfluid(2,c,b,a) = v_2(a,b,c) * factor
          vfluid(3,c,b,a) = v_3(a,b,c) * factor

          if (enable_gr) then
            phiconf(c,b,a) = phi(a,b,c)
            alpha  (c,b,a) = alpha_in(a,b,c)
            beta (1,c,b,a) = beta_1(a,b,c)
            beta (2,c,b,a) = beta_2(a,b,c)
            beta (3,c,b,a) = beta_3(a,b,c)
          else
            phiconf(c,b,a) = 1.0
            alpha  (c,b,a) = 1.0
            beta (:,c,b,a) = 0.0
          endif
        enddo
      enddo
    enddo

    ! Get f_eq, kappa_a, and kappa_s
    call calculate_opacities(rho, t, xnnu, nx)

    ! Determine jump locations
    call detect_jumps(rho)

    ! Boost all cells for change in fluid velocity and metric
    call time_boost

    ! Reset hydro source terms
    qye = 0.0
    qen = 0.0
    qmo = 0.0

    ! Perform substeps
    do istep = 1, nsteps
      call step(dt)

      ! Calculate source terms in comoving frame
      call hydro_backreaction(qye, qen, qmo)
    enddo

    ! Convert units of the source terms
    ! Sign flip because this is the source term for the fluid
    ! Average the source terms over the substeps
    const = 1.0 / (2.0*pc_pi*wc_hc)**3
    constc = const * pc_cl ! factor of c from the opacity term
    insteps = 1.0 / real(nsteps)
    qye    = - qye * pc_mb     * constc * insteps
    qen    = - qen * pc_meverg * constc * insteps
    qmo    = - qmo * pc_meverg * const  * insteps

    ! Copy the momentum components into the output arrays
    qmo_r = qmo(1,:,:,:)

    call additional_terms(dnu, enu, fnu, pnu, enu_lab, fnu_lab)

    ! Add units
    dnu = dnu * const
    enu = enu * const  * pc_meverg
    fnu = fnu * constc * pc_meverg
    pnu = pnu * const  * pc_meverg

    if (metric_contribution) then
      ! Geometric units for neutrino contribution to metric
      enu_lab = enu_lab * const * pc_meverg * pc_gc / pc_cl**4
      fnu_lab = fnu_lab * const * pc_meverg * pc_gc / pc_cl**4

      ! Set central fnu = 0
      fnu_lab(1,:,:) = 0.0

    else
      enu_lab = 0.0
      fnu_lab = 0.0
    endif

  end subroutine boltzmann_step

  subroutine write_boltzmann_output(filename, nstep_in)
    use rice_config, only: use_boltzmann
    use rice_output, only: write_output
    use rice_grid,   only: nstep
    character(200), intent(in) :: filename
    integer,        intent(in) :: nstep_in

    nstep = nstep_in

    ! Only write output if BT solver is being used
    if (use_boltzmann) call write_output(filename)

  end subroutine write_boltzmann_output

  subroutine read_boltzmann_output(filename)
    use rice_config, only: use_boltzmann
    use rice_output, only: read_output

    character(200), intent(in) :: filename

    ! Only read output if BT solver is being used
    if (use_boltzmann) call read_output(filename)

  end subroutine read_boltzmann_output

  subroutine calculate_opacities(rho, t, xnnu, nx)
    use rice_config, only: neps
    use rice_grid,   only: eps, f_eq, kappa_a, kappa_s, &
                           a_s, a_e, b_s, b_e, c_s, c_e
    use eos_sn2,     only: eos

    ! CoCoNuT indexing is the other way around
    integer, intent(in)   :: nx ! number of species
    real,    intent(in)   :: rho (a_s:a_e,b_s:b_e,c_s:c_e)  ! density in geometric units
    real,    intent(in)   :: t   (a_s:a_e,b_s:b_e,c_s:c_e)  ! temperature in Kelvin
    real,    intent(in)   :: xnnu(a_s:a_e,b_s:b_e,c_s:c_e,1:nx) ! mass fraction of species

    real :: den  (a_s:a_e)
    real :: nby  (a_s:a_e)
    real :: tem  (a_s:a_e)
    real :: tmev (a_s:a_e)
    real :: tmev1(a_s:a_e)
    real :: xiq  (a_s:a_e,1:nx)
    real :: etanp(a_s:a_e)
    real :: etapn(a_s:a_e)
    real :: etann(a_s:a_e)
    real :: etapp(a_s:a_e)
    real :: dup  (a_s:a_e)
    real :: dnuc (a_s:a_e)

    real :: cuq(a_s:a_e)
    real :: ceq(a_s:a_e)
    real :: cnq(a_s:a_e)
    real :: cpq(a_s:a_e)

    ! For the EOS call
    real :: xhrepq(a_s:a_e)
    real :: zah   (a_s:a_e,2)
    real :: ein   (a_s:a_e)
    real :: pre   (a_s:a_e)
    real :: gamc  (a_s:a_e)
    real :: s     (a_s:a_e)
    real :: eos_self(2)
    real :: eos_children(2)
    logical :: ler

    integer :: abrems

    integer:: a, b, c, i

    ! Over rays (array ordering is sub-optimal for now)
    do c = c_s, c_e
      do b = b_s, b_e
        ! Convert density to grams
        den(:) = rho(:,b,c) * pc_geog

        ! Get number of baryons
        nby(:) = den(:) / pc_mb

        ! Convert temperature to MeV
        tem(:)   = t(:,b,c)
        tmev(:)  = tem(:) * pc_kmev
        tmev1(:) = 1.0 / tmev(:)

        ! Repack species array
        xiq(:,:) = xnnu(:,b,c,:)

        etanp(:) = nby(:) * xiq(:,1)
        etapn(:) = nby(:) * xiq(:,2)
        etann(:) = nby(:) * xiq(:,1)
        etapp(:) = nby(:) * xiq(:,2)

        ! Call equation of state to get chemical potentials
        call eos(den(:), tem(:), xiq(:,:), &
                 xhrepq(:), zah(:,:), ein(:), pre(:), gamc(:), s(:), & ! unused variables
                 cuq(:), ceq(:), cnq(:), cpq(:), & ! chemical potentials
                 eos_self(:), eos_children(:), mode=1, nsemode=1, ler=ler)

        dup(:) = 0.0
        ! Apply blocking correction
        do a = a_s, a_e
          call blocking_correction(tmev(a), cnq(a), cpq(a), &
                                   etanp(a), etapn(a), etann(a), etapp(a), dup(a))
        enddo

        where (zah(:,2) /= 0.0)
          dnuc(:) = den(:) * xhrepq(:) / (pc_mb * zah(:,2))
        else where
          dnuc(:) = 0.0
        end where

        ! Determine where we need to calculate Bremsstrahlung
        abrems = 1
        ! This still works for radial decomposition because
        ! abrems will get pushed out to the end
        do a = min(a_s, 2), a_e
          if (den(a) > denbrems) then
            abrems = a
          endif
        enddo

        ! Over each energy bin
        do i = 1, neps
          ! Over spatial coordinate
          do a = a_s, a_e
            ! Calculate scattering coefficient
            kappa_s(:,i,c,b,a) = kappa_s_all(eps(i), tmev(a), etapp(a), etann(a), nby(a), &
                                 nx, xiq(a,:), zah(a,:), dnuc(a))

            ! Calculate equilibrium values
            f_eq(1,i,c,b,a) = f_eq_en (eps(i), cuq(a), tmev1(a))
            f_eq(2,i,c,b,a) = f_eq_ean(eps(i), cuq(a), tmev1(a))
            f_eq(3,i,c,b,a) = f_eq_mt (eps(i),         tmev1(a))

            ! Calculate absorption coefficient
            kappa_a(1,i,c,b,a) = kappa_a_en (eps(i), tmev1(a), cuq(a), ceq(a), etanp(a), dup(a), &
                                             zah(a,:), dnuc(a))
            kappa_a(2,i,c,b,a) = kappa_a_ean(eps(i), tmev1(a), cuq(a), ceq(a), etapn(a), dup(a))
            kappa_a(3,i,c,b,a) = kappa_a_mt(eps(i), tmev(a), tmev1(a), etapp(a), etann(a), &
                                            kappa_s(3,i,c,b,a), (a <= abrems))

          enddo
        enddo

      enddo
    enddo

  end subroutine calculate_opacities

  function f_eq_en(eps, cuq, tmev1) result(f_eq)
    real, intent(in) :: eps
    real, intent(in) :: cuq
    real, intent(in) :: tmev1

    real :: f_eq

    f_eq = bnorm / (1.0 + safe_exp((eps - cuq)*tmev1))

  end function f_eq_en

  function f_eq_ean(eps, cuq, tmev1) result(f_eq)
    real, intent(in) :: eps
    real, intent(in) :: cuq
    real, intent(in) :: tmev1

    real :: f_eq

    ! Sign difference                  v
    f_eq = bnorm / (1.0 + safe_exp((eps + cuq)*tmev1))
  end function f_eq_ean

  function f_eq_mt(eps, tmev1) result(f_eq)
    real, intent(in) :: eps
    real, intent(in) :: tmev1

    real :: f_eq

    f_eq = bnorm / (1.0 + safe_exp(eps*tmev1))
  end function f_eq_mt

  function kappa_a_en(eps, tmev1, cuq, ceq, etanp, dup, zah, dnuc) result(kappa_a)
    real, intent(in) :: eps
    real, intent(in) :: tmev1
    real, intent(in) :: cuq
    real, intent(in) :: ceq
    real, intent(in) :: etanp
    real, intent(in) :: dup
    real, intent(in) :: zah(2)
    real, intent(in) :: dnuc

    real, parameter :: delta = 3.0

    real :: eargn, earge
    real :: fn, fe
    real :: scre
    real :: con_en

    real :: fnp, fnh
    real :: qprim
    real :: scr2, scr3

    real :: kappa_a_neutrons = 0.0
    real :: kappa_a_nuclei = 0.0
    real :: kappa_a

    ! Absorption on neutrons
    eargn = (eps - cuq)*tmev1
    earge = (eps - ceq + wc_mq + dup)*tmev1
    fn = 1.0 + safe_exp(-eargn)
    if (earge > 35.0) then
      fe = 1.0
    else
      fe = safe_exp(earge) / (1.0 + safe_exp(earge))
    endif
    scre = 0.25*wc_s0*max(eps + wc_mq + dup, 0.0) &
           / wc_me**2 * sqrt(max(max(eps + wc_mq + dup, 0.0)**2 - wc_me**2, 0.0))
    con_en = (1.0 + 3.0*wc_ga2)*scre
    kappa_a_neutrons = con_en*fe*fn*etanp

    ! Absorption on nuclei
    fnp = max(zah(1) - 20.0, 0.0)
    fnp = min(fnp, 8.0)
    fnh = max(40.0 - (zah(2) - zah(1)), 0.0)
    fnh = min(fnh, 6.0)

    scr3 = fnp * fnh

    ! Numerator is different for nuclei
    ! For absorption onto neutrons, Q  = m_n *c^2 - m_p *c^2
    ! For absorption onto nuclei,   Q' = mu_n*c^2 - mu_p*c^2 + delta
    qprim = ceq - cuq + delta
    earge = (eps + qprim - ceq)*tmev1

    if (eargn > 35.0 .and. earge > 35.0) then
       scr2 = safe_exp(eargn - earge)
    else
       scr2 = (1.0 + safe_exp(eargn)) / (1.0 + safe_exp(earge))
    endif
    if (eps + qprim <= wc_me) then
       scr2 = 0.0
    else
       scr2 = scr2 * wc_s0*wc_ga2/14.0 * (eps + qprim) / wc_me**2 * sqrt((eps + qprim)**2 - wc_me**2)
    endif
    kappa_a_nuclei = scr2 * scr3 * dnuc

    ! kappa_a = pc_cl * (kappa_a_neutrons + kappa_a_nuclei)
    kappa_a = kappa_a_neutrons + kappa_a_nuclei

  end function kappa_a_en

  function kappa_a_ean(eps, tmev1, cuq, ceq, etapn, dup) result(kappa_a)
    real, intent(in) :: eps
    real, intent(in) :: tmev1
    real, intent(in) :: cuq
    real, intent(in) :: ceq
    real, intent(in) :: etapn
    real, intent(in) :: dup

    real :: eargn, earge
    real :: fn, fe
    real :: scre
    real :: con_en

    real :: kappa_a

    ! Absorption on protons
    ! This sign is different
    !            v
    eargn = (eps + cuq)*tmev1
    earge = (eps + ceq - wc_mq - dup)*tmev1
    fn = 1.0 + safe_exp(-eargn)
    if (earge > 35.0) then
      fe = 1.0
    else
      fe = safe_exp(earge) / (1.0 + safe_exp(earge))
    endif
    ! These signs are different
    !                         v
    scre = 0.25*wc_s0*max(eps - wc_mq - dup, 0.0) &
    / & !                         v
    wc_me**2 * sqrt(max(max(eps - wc_mq - dup, 0.0)**2 - wc_me**2, 0.0))
    con_en = (1.0 + 3.0*wc_ga2)*scre

    ! kappa_a = pc_cl*con_en*fe*fn*etapn
    kappa_a = con_en*fe*fn*etapn
  end function kappa_a_ean

  function kappa_a_mt(epsi, tmev, tmev1, etapp, etann, kappa_s, brems) result(kappa_a)
    use rice_config, only: neps
    use rice_grid,   only: eps, eps_if

    real,    intent(in) :: epsi
    real,    intent(in) :: tmev
    real,    intent(in) :: tmev1
    real,    intent(in) :: etapp
    real,    intent(in) :: etann
    real,    intent(in) :: kappa_s
    logical, intent(in) :: brems

    real :: scr2, scr3
    real :: cbrems
    real :: cexc
    real :: kappa_a

    integer :: i

    ! Bremsstrahlung
    cbrems = 0.0
    if (brems) then
      do i = 1, neps
        ! alpha does not appear in these terms (unlike in FMT)
        ! because the energy grid is same in every cell
        scr2 = (epsi + eps(i)) * tmev1
        scr3 = 1.0 / (1.0 + safe_exp(eps(i) * tmev1))

        cbrems = cbrems + (eps_if(i)**3 - eps_if(i-1)**3) / 3.0 &
                / (epsi + eps(i)) * sqrt(0.5 * pc_pi / scr2) &
                * ((1.0 - scr3) * safe_exp(-scr2) + scr3)

      enddo

      ! cbrems = pc_cl * cbrems * 0.25 * 0.5 * wc_s0 * wc_ga2 &
      ! / sqrt(pc_pi)**5 * (1.0 / 135.0)**4 * (wc_mn/wc_me)**2 &
      ! / sqrt(wc_mn * tmev) * wc_hc**3 * (etapp**2 + etann**2 + 28.0/3.0*etapp*etann)
      cbrems = cbrems * 0.25 * 0.5 * wc_s0 * wc_ga2 &
      / sqrt(pc_pi)**5 * (1.0 / 135.0)**4 * (wc_mn/wc_me)**2 &
      / sqrt(wc_mn * tmev) * wc_hc**3 * (etapp**2 + etann**2 + 28.0/3.0*etapp*etann)
    endif

    ! Energy exchange in scattering reactions
    cexc = kappa_s * abs(epsi - 3.0*tmev) / wc_mnc

    kappa_a = cbrems + cexc

  end function kappa_a_mt

  function kappa_s_all(epsi, tmev, etapp, etann, nby, nx, xnnu, zah, dnuc) result(kappa_s)
    use nucparam, only: pc_nuc
    real,    intent(in) :: epsi
    real,    intent(in) :: tmev
    real,    intent(in) :: etapp
    real,    intent(in) :: etann
    real,    intent(in) :: nby
    integer, intent(in) :: nx
    real,    intent(in) :: xnnu(1:nx)
    real,    intent(in) :: zah(2)
    real,    intent(in) :: dnuc

    real, parameter :: ga_s = -0.1 ! nucleon strangeness

    real :: scre
    real :: scr2
    real :: ye
    real :: s_a

    real :: kappa_s_nucleons
    real :: kappa_s_nuclei
    real :: kappa_s

    integer :: i_n

    ! Horowitz (2016) spin response
    scr2 = nby * 1.e-39
    ye = xnnu(nx)
    s_a = 1.0 / (1.0 +                        &
          920.0 * scr2 *                      & !A
          (1.0 - ye + ye**2) /                &
          tmev**1.22 *                        &
          (1.0 + 3.05 / sqrt(sqrt(tmev))**3 / & !B
          safe_exp(6140.0 * scr2 * ye *       & !C
          (1.0 - ye) / sqrt(tmev) +           &
          1.5e13 * scr2**4 / tmev** 6)))

    ! Scattering on nucleons
    scre = 0.25 * wc_s0 * (epsi/wc_me)**2

    ! No weak magnetism
    kappa_s_nucleons = scre * (                                                      &
                etapp * ((1.0 - wc_cv)**2 + 0.75 * ( sqrt(wc_ga2) - ga_s)**2  * s_a) &
              + etann * (      0.25 * (1.0 + 3.0 * (-sqrt(wc_ga2) - ga_s)**2) * s_a) &
              )

    ! Scattering on nuclei
    kappa_s_nuclei = 0.0
    ! if (zah(2) == 0.0) then ! low density regime
    !
    !   do i_n = 3, nx-1
    !     scre = 0.125*wc_s0*(epsi/wc_me)**2
    !     scr2 = nby * xnnu(i_n) / pc_nuc(i_n,2)
    !     scre = scre * 0.25 * pc_nuc (i_n,2)**2 *        &
    !            (wc_ca - wc_cv + (2.0 - wc_ca - wc_cv) * &
    !            (2.0 * pc_nuc(i_n,1) - pc_nuc(i_n,2)) /  &
    !            pc_nuc(i_n,2))**2
    !     kappa_s_nuclei = kappa_s_nuclei + scr2 * scre
    !   enddo
    !
    ! else ! high density regime
    !
    !   scre = 0.125*wc_s0*(epsi/wc_me)**2
    !   scr2 = dnuc
    !   kappa_s_nuclei = kappa_s_nuclei + scr2 * scre *           &
    !                    0.25 * zah(2)**2 *                       &
    !                    (wc_ca - wc_cv + (2.0 - wc_ca - wc_cv) * &
    !                    (2.0 * zah(1) - zah(2)) / zah(2))**2
    !
    ! endif

    ! kappa_s = pc_cl * (kappa_s_nucleons + kappa_s_nuclei)
    kappa_s = kappa_s_nucleons + kappa_s_nuclei

  end function kappa_s_all

  function safe_exp(x) result(expx)
    use machcons, only: bigeexp
    real, intent(in)  :: x
    real :: expx

    expx = exp(min(x, bigeexp))
  end function safe_exp

  subroutine blocking_correction(tmev, cnq, cpq, etanp, etapn, etann, etapp, dup)
    use mod_degeneracy_parameter, only: etapike
    use machcons,                 only: bigeexp
    real, intent(in)    :: tmev
    real, intent(in)    :: cnq
    real, intent(in)    :: cpq
    real, intent(inout) :: etanp
    real, intent(inout) :: etapn
    real, intent(inout) :: etann
    real, intent(inout) :: etapp
    real, intent(inout) :: dup

    real, parameter :: fac = 3.0*pc_pi**2*wc_hc**3
    real, parameter :: pcon = 0.5*wc_hc*wc_hc/wc_mb
    real, parameter :: pi32 = 3.0*pc_pi*pc_pi

    real :: scr2, scr3, escr
    real :: etnucn, etnucp, detnucnp

    ! etapn and etanp
    scr2 = fac / sqrt(2.0*wc_mnc*tmev)**3 * max(etanp, 1.0e-16)
    scr3 = fac / sqrt(2.0*wc_mnc*tmev)**3 * max(etapn, 1.0e-16)

    etnucn = min(etapike(scr2), bigeexp)
    etnucp = min(etapike(scr3), bigeexp)
    detnucnp = etnucn - etnucp

    dup = (cnq - cpq) - tmev*detnucnp

    if (abs(detnucnp) > 1.0e-7) then
      escr = exp(detnucnp)
      scr2 = escr - 1.0
      scr3 = 1.0/escr - 1.0
    else
      scr2 =  detnucnp + 0.5*detnucnp**2
      scr3 = -detnucnp + 0.5*detnucnp**2
    endif

    scr2 = sign(max(abs(scr2),1.0e-15), scr2)
    scr3 = sign(max(abs(scr3),1.0e-15), scr3)

    scr2 = (etanp - etapn) / scr2
    scr3 = (etapn - etanp) / scr3

    etapn = max(min(etapn, scr2), 0.0)
    etanp = max(min(etanp, scr3), 0.0)

    ! etann and etapp
    scr2 = max(etapp, 1.0e-15) ! 1e-15 protons/ccm is a small number
    scr2 = 1.5*tmev / (pcon*(scr2 * pi32)**(2.0/3.0))
    etapp = etapp * scr2 / sqrt(1.0 + scr2**2)

    scr3 = max(etann, 1.0e-15) ! 1e-15 neutrons/ccm is a small number
    scr3 = 1.5*tmev / (pcon*(scr3 * pi32)**(2.0/3.0))
    etann = etann * scr3 / sqrt(1.0 + scr3**2)

  end subroutine blocking_correction

  subroutine hydro_backreaction(qye, qen, qmo)
    use rice_config, only: nflav, npsi, nmu, neps
    use rice_grid,   only: f, f_eq, kappa_a, kappa_s, &
                           a_s, a_e, b_s, b_e, c_s, c_e
    real, intent(inout) :: qye(    a_s:a_e,b_s:b_e,c_s:c_e)
    real, intent(inout) :: qen(    a_s:a_e,b_s:b_e,c_s:c_e)
    real, intent(inout) :: qmo(1:3,a_s:a_e,b_s:b_e,c_s:c_e)

    integer :: a, b, c

    !$omp parallel do
    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e
          call source_terms(f      (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a), &
                            f_eq   (1:nflav,1:neps,                c,b,a), &
                            kappa_a(1:nflav,1:neps,                c,b,a), &
                            kappa_s(1:nflav,1:neps,                c,b,a), &
                            qye    (    a,b,c)      , &
                            qen    (    a,b,c)      , &
                            qmo    (1:3,a,b,c)      )
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine hydro_backreaction

  subroutine source_terms(f, f_eq, kappa_a, kappa_s, qye_tot, qen_tot, qmo_tot)
    use rice_config, only: nflav, npsi, nmu, neps
    use rice_grid,   only: domega, domega_scaled, eps, eps_if, u_com, epsvol, mu, &
                           eps2domegadeps, eps3domegadeps
    real, intent(in)    :: f       (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(in)    :: f_eq    (1:nflav,1:neps                )
    real, intent(in)    :: kappa_a (1:nflav,1:neps                )
    real, intent(in)    :: kappa_s (1:nflav,1:neps                )
    real, intent(inout) :: qye_tot
    real, intent(inout) :: qen_tot
    real, intent(inout) :: qmo_tot(3)

    real :: qye, qye_en, qye_ean
    real :: qen, qen_en, qen_ean, qen_mt
    real :: qmo_en(3), qmo_ean(3), qmo_mt(3)

    real :: ks_int(1:nflav,1:neps)

    real :: c
    real :: cd
    real :: epscd

    real :: const, constc

    integer :: i, j, k, l

    ! Calculate term for scattering integral
    ks_int(:,:) = 0.0
    do i = 1, neps
      do l = 1, nflav
        do j = -nmu, nmu
          do k = 1, npsi
            ks_int(l,i) = ks_int(l,i) + f(l,i,k,j) * domega_scaled(k,j)
          enddo
        enddo
      enddo
    enddo

    qye_en  = 0.0
    qye_ean = 0.0

    qen_en  = 0.0
    qen_ean = 0.0
    qen_mt  = 0.0

    qmo_en  = 0.0
    qmo_ean = 0.0
    qmo_mt  = 0.0

    do i = 1, neps
      do l = 1, nflav
        do j = -nmu, nmu
          do k = 1, npsi

            ! Collision integral
            c = kappa_a(l,i) * (f_eq(l,i)   - f(l,i,k,j)) &
              + kappa_s(l,i) * (ks_int(l,i) - f(l,i,k,j))

            ! l==3 represents 4 flavours of neutrinos --> contributes 4 times
            if (l == 3) then
              c = c * 4.0
            endif

            cd    = c*eps2domegadeps(i,k,j)
            epscd = c*eps3domegadeps(i,k,j)

            ! Lepton number source term
            if (l == 1) then
              qye_en    = qye_en    + cd
              qen_en    = qen_en    + epscd
              qmo_en(:) = qmo_en(:) + u_com(1:3,k,j)*epscd
            elseif (l == 2) then
              qye_ean    = qye_ean    - cd
              qen_ean    = qen_ean    + epscd
              qmo_ean(:) = qmo_ean(:) + u_com(1:3,k,j)*epscd
            elseif (l == 3) then
              qen_mt    = qen_mt    + epscd
              qmo_mt(:) = qmo_mt(:) + u_com(1:3,k,j)*epscd
            endif

          enddo
        enddo
      enddo
    enddo

    qye_tot    = qye_tot + qye_en + qye_ean
    qen_tot    = qen_tot + qen_en + qen_ean + qen_mt
    qmo_tot    = qmo_tot + qmo_en + qmo_ean + qmo_mt

  end subroutine source_terms

  subroutine additional_terms(dnu, enu, fnu, pnu, enu_lab, fnu_lab)
    use rice_config, only: nflav, npsi, nmu, neps
    use rice_grid,   only: f, f_eq, kappa_a, kappa_s, domega, eps, epsvol, mu, &
                           a_s, a_e, b_s, b_e, c_s, c_e, &
                           phiconf, vfluid, &
                           eps2domegadeps, eps3domegadeps, &
                           mueps3domegadeps, mu2eps3domegadeps
    use rice_sr,     only: lambda_transform

    real, intent(out) :: dnu    (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real, intent(out) :: enu    (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real, intent(out) :: fnu    (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real, intent(out) :: pnu    (a_s:a_e,b_s:b_e,c_s:c_e,1:nflav)
    real, intent(out) :: enu_lab(a_s:a_e,b_s:b_e,c_s:c_e) ! lab only needs total
    real, intent(out) :: fnu_lab(a_s:a_e,b_s:b_e,c_s:c_e)

    integer :: a, b, c, i, j, k, l

    real :: lambda(0:3,0:3)
    real :: lambda_t(0:3)
    real :: t(0:3,0:3)

    dnu = 0.0
    ! enu = 0.0
    ! fnu = 0.0
    ! pnu = 0.0

    enu_lab = 0.0
    fnu_lab = 0.0

    !$omp parallel &
    !$omp default(none) &
    !$omp private(a,b,c,l,j,k,i) &
    !$omp private(t, lambda, lambda_t) &
    !$omp shared(a_s, a_e, b_s, b_e, c_s, c_e) &
    !$omp shared(f, vfluid, phiconf) &
    !$omp shared(eps2domegadeps, eps3domegadeps, mueps3domegadeps, mu2eps3domegadeps) &
    !$omp shared(dnu, enu, fnu, pnu) &
    !$omp shared(enu_lab, fnu_lab)
    !$omp do
    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e

          do l = 1, nflav

            do j = -nmu, nmu
              do k = 1, npsi
                do i = 1, neps
                  dnu(a,b,c,l) = dnu(a,b,c,l) + f(l,i,k,j,c,b,a) * eps2domegadeps(i,k,j)

                  ! The following can be computed using the energy-stress tensor
                  ! enu(a,b,c,l) = enu(a,b,c,l) + f(l,i,k,j,c,b,a) * eps3domegadeps(i,k,j)
                  ! fnu(a,b,c,l) = fnu(a,b,c,l) + f(l,i,k,j,c,b,a) * mueps3domegadeps(i,k,j)
                  ! pnu(a,b,c,l) = pnu(a,b,c,l) + f(l,i,k,j,c,b,a) * mu2eps3domegadeps(i,k,j)
                enddo
              enddo
            enddo

            t = stress_energy_tensor(a,b,c,l)

            ! Comoving frame
            enu(a,b,c,l) = t(0,0)
            fnu(a,b,c,l) = t(1,0)
            pnu(a,b,c,l) = (t(1,1) + t(2,2) + t(3,3)) * third

            ! Compute lab frame energy density and flux
            lambda = lambda_transform(vfluid(:,c,b,a))
            lambda_t = matmul(t, lambda(:,0))
            enu_lab(a,b,c) = enu_lab(a,b,c) + phiconf(c,b,a)**6 * sum(lambda(:,0)*lambda_t)
            fnu_lab(a,b,c) = fnu_lab(a,b,c) + phiconf(c,b,a)**8 * sum(lambda(:,1)*lambda_t)

          enddo

        enddo
      enddo
    enddo
    !$omp enddo
    !$omp end parallel

  end subroutine additional_terms

  function stress_energy_tensor(a,b,c,l) result(t)
    use rice_config, only: nmu, npsi, neps, nflav
    use rice_grid,   only: f, eps3domegadeps, u_com

    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    integer, intent(in) :: l

    real :: t(0:3,0:3)
    real :: fe

    integer :: i, j, k

    t = 0.0

    do j = -nmu, nmu
      do k = 1, npsi
        do i = 1, neps
          fe = f(l,i,k,j,c,b,a) * eps3domegadeps(i,k,j)

          t(0,0) = t(0,0) + fe
          t(1,0) = t(1,0) + fe*u_com(1,k,j)
          t(2,0) = t(2,0) + fe*u_com(2,k,j)
          t(3,0) = t(3,0) + fe*u_com(3,k,j)

          t(0,1) = t(0,1) + fe*u_com(1,k,j)
          t(1,1) = t(1,1) + fe*u_com(1,k,j)*u_com(1,k,j)
          t(2,1) = t(2,1) + fe*u_com(2,k,j)*u_com(1,k,j)
          t(3,1) = t(3,1) + fe*u_com(3,k,j)*u_com(1,k,j)

          t(0,2) = t(0,2) + fe*u_com(2,k,j)
          t(1,2) = t(1,2) + fe*u_com(1,k,j)*u_com(2,k,j)
          t(2,2) = t(2,2) + fe*u_com(2,k,j)*u_com(2,k,j)
          t(3,2) = t(3,2) + fe*u_com(3,k,j)*u_com(2,k,j)

          t(0,3) = t(0,3) + fe*u_com(3,k,j)
          t(1,3) = t(1,3) + fe*u_com(1,k,j)*u_com(3,k,j)
          t(2,3) = t(2,3) + fe*u_com(2,k,j)*u_com(3,k,j)
          t(3,3) = t(3,3) + fe*u_com(3,k,j)*u_com(3,k,j)
        enddo
      enddo
    enddo

  end function stress_energy_tensor

  subroutine detect_jumps(rho)
    use rice_grid,   only: jump, a_s, a_e, b_s, b_e, c_s, c_e
    real, intent(in) :: rho(a_s:a_e,b_s:b_e,c_s:c_e) ! in geometric units

    real :: kappa_l_min, kappa_l_max, kappa_r_min, kappa_r_max

    integer :: a, b, c

    jump = .false.

    do c = c_s, c_e
      do b = b_s, b_e
        do a = a_s, a_e-1
          if (rho(a,b,c) > 2.0*rho(a+1,b,c) .or. rho(a,b,c) < 0.5*rho(a+1,b,c)) then
            ! Set jump based on jump in r direction, even for lateral interfaces
            jump(c-1,b-1,a) = .true.
            jump(c-1,b  ,a) = .true.
            jump(c,  b-1,a) = .true.
            jump(c,  b  ,a) = .true.
          endif
        enddo
      enddo
    enddo

  end subroutine detect_jumps

end module rice_coconut_interface
