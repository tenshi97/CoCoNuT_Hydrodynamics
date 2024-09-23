program boltzmann
  use rice_config, only: tmax, noutputs, stepmax, heun
  use rice_tune,   only: tune_areas
  use rice_step,   only: step
  use rice_dt,     only: get_dt
  use rice_ic,     only: setup_ic
  use rice_output, only: write_output, init_hdf5, finish_hdf5
  use rice_timing, only: start_timing, stop_timing, print_timing
  use rice_grid,   only: setup_grid, time, index, nstep, f, f_old
  use rice_mpi,    only: setup_mpi, finish_mpi, calculate_ncyc

  implicit none
  !  -1 0  1   ...  n  n+1
  !  |  |  |        |  |
  !   0  1  2 ... n  n+1

  real :: dt, dtoutput, nextoutput

  call setup_mpi(init=.true.)
  call init_hdf5
  call setup_grid
  call calculate_ncyc

  index = -1
  call setup_ic
  if (index == -1) call write_output

  ! dt requirement doesn't change
  call get_dt(dt)

  ! Determine output dt
  dtoutput = tmax / real(noutputs)
  nextoutput = time + dtoutput

  if (heun) then
    print*, "(main) Using Heun's method"
  else
    print*, "(main) Using first order time integration"
  endif

  do while (time < tmax .and. nstep < stepmax)

    ! Time the first few steps
    if (nstep == 1) call start_timing
    if (nstep == 3) then
      call stop_timing
      call print_timing(3, int(tmax/dt))
    endif

    if (heun) then ! second order
      f_old = f
      call step(dt)
      call step(dt)
      f = 0.5 * (f_old + f)

    else ! first order
      call step(dt)
    endif

    nstep = nstep + 1
    if (time >= nextoutput) then
        call write_output
        nextoutput = nextoutput + dtoutput
    endif

    time = time + dt

  enddo

  ! Write final output
  call write_output
  call finish_hdf5

  call finish_mpi

  print*,'(main) Finished!'

end program boltzmann
