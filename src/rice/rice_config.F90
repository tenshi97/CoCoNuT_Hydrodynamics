! WARNING: This file will be overwritten by the test suite

module rice_config

    use rice_constants, only: pi
    implicit none

    ! Grid size
    integer :: nr = 100                     ! number of radial cells
    integer :: ntheta = 1                   ! number of altitudinal cells
    integer :: nphi = 1                     ! number of azimuthal cells
    integer, parameter :: neps = 12          ! number of cells on the energy grid

    logical :: uniform_mu = .false.         ! space the momentum space grid uniformly between -1 and 1 instead of using the Gauss-Lobatto colocation points
    integer, parameter :: nmu = 4           ! number of cells in the mu direction on the momentum grid
    integer, parameter :: npsi = 4          ! number of cells in the angular direction on the momentum grid
    integer, parameter :: nflav = 3         ! number of particle flavours

    ! Grid dimensions
    logical :: auto_rmin = .false.          ! automatically set rmin to rmax**(1/nr) - 1
    real :: rmin = 0.1                      ! radial coordinate of second-innermost radial interface (the innermost is at zero)
    real :: rmax = 6.0                      ! radial coordinate of the outermost radial interface

    ! Automatically set dtheta and dphi, ignore values below
    logical :: auto_grid = .true.           ! automatically set theta and phi interfaces for cells of reasonable aspect ratio

    real :: thetamin = 1.5053464798         ! lower bound of theta wedge
    real :: thetamax = 1.6362461737         ! upper bound of theta wedge
    real :: phimin = 0.0                    ! lower bound of phi wedge
    real :: phimax = 0.1308996939           ! upper bound of phi wedge
    real :: epsmin = 0.0                    ! lower bound of energy grid
    real :: epsmax = 1.0                    ! upper bound of energy grid

    ! Check grid aspect ratio
    logical :: check_grid = .false.         ! check that the aspect ratio of cells is sensible

    ! Time integration
    real :: tmax = 10.0                     ! physical end time
    real :: cfl_factor = 0.8                ! CFL factor calculated using smallest cell
    real :: stepmax = 1000000               ! maximum number of time steps

    ! Option to disable fluxes in each direction
    logical :: flux_r = .true.              ! enable radial flux
    logical :: flux_theta = .true.          ! enable theta flux
    logical :: flux_phi = .true.            ! enable phi flux

    ! Second order time integration
    logical :: heun = .false.               ! use heun's method for time integration

    ! Do rotation in transform
    logical :: transform_rotate = .true.    ! account for rotation between cells during transformation

    ! Perform geometric rotation for Lax-Wendroff
    logical :: lw_rotate = .true.           ! account for rotation between cells during LW reconstruction

    ! LW flux only
    logical :: lw_only = .false.            ! only use LW fluxes

    ! Reconstruct f*r^2 along the r direction
    logical :: reconstruct_fr2 = .false.    ! perform reconstruction using f*r^2 along the radial direction

    ! Slope limiter (0: no reconstruction, 1: minmod, 2: MC)
    integer :: slope_limiter = 2            ! slope limiter for reconstruction

    ! Mu advection everywhere
    logical :: mu_advection = .false.       ! advect the distribution function in mu everywhere on each timestep

    ! Core treatments
    integer :: nr_nolateral = 0             ! number of central zones to disable lateral flux
    integer :: nr_isotropic = 0             ! number of central zones to make isotropic
    integer :: nr_diffusive = 0             ! number of central zones to use diffusion approximation
    integer :: nr_evolve_mu = 0             ! number of central zones to evolve mu (advect according to timestep)
    logical :: forward_inner = .false.      ! convert innermost zone to mu>0
    logical :: scatter_inner = .false.      ! redistribute mu=-1 in innermost zone to adjacent mu
    logical :: average_inner = .false.      ! average +mu and -mu at core
    logical :: average_poles = .false.      ! average psi at poles
    logical :: square_inner = .false.       ! make innermost interface area non-zero
    logical :: isotropic_3d = .false.       ! make central zone isotropic and average across other central zones
    logical :: cartoon_grid = .false.       ! use cartoon grid in 2d
    logical :: cartoon_rotate = .false.     ! rotate BCs for cartoon grid

    ! Physics
    real :: clight = 1.0                    ! speed of light

    ! Outputs
    integer :: noutputs = 10                ! number of outputs during the entire run

    ! Boundary conditions
    logical :: closed_boundaries = .false.  ! disable outflow from the grid

    ! Testing
    logical :: ic_reset = .false.           ! reset index, nstep, and time to zero in output (for generating ICs)

    ! MPI domain decomposition
    integer :: mpi_nr = 1                   ! size of MPI grid in radial direction
    integer :: mpi_ntheta = 1               ! size of MPI grid in theta direction
    integer :: mpi_nphi = 1                 ! size of MPI grid in phi direction

    ! Give identical results regardless of number of MPI tasks
    logical :: mpi_identical = .false.      ! Runs slower, for debugging and testing only
    logical :: zero_beta = .false.          ! Set the shift vector to zero

    ! CoCoNuT coupling
#ifdef CFC_TRANSPORT2
    logical :: use_boltzmann = .false.
#else
    logical :: use_boltzmann = .false.       ! use the BT solver (uses the FMT solver if disabled)
#endif
    logical :: metric_contribution = .true. ! neutrino contribution to the metric
    logical :: enable_gr = .true.           ! use GR for radiation, otherwise set alpha=phi=1 and beta=0



end module rice_config
