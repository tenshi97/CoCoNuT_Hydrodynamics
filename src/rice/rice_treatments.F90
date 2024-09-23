module rice_treatments
  implicit none

contains

  subroutine evolve_mu(r, dt, u_mu, kappa_a, kappa_s)
    use rice_config,    only: neps, nmu, npsi, nflav, clight, mu_advection
    use rice_grid,      only: mu, u_com
    real, intent(in)    :: r
    real, intent(in)    :: dt
    real, intent(out)   :: u_mu(0:3,1:npsi,-nmu:nmu)
    real, optional, intent(in) :: kappa_a(1:nflav,1:neps)
    real, optional, intent(in) :: kappa_s(1:nflav,1:neps)

    real :: cdt
    real :: r_new
    real :: r2cdt2
    real :: mu_new
    real :: rescale

    integer :: j

    if (.not. mu_advection) then
      u_mu = u_com
      return
    endif

    ! If not using a real dt, evolve to prevent central buildup
    if (present(kappa_a) .and. present(kappa_s)) then
      cdt = r * exp(-clight*dt*maxval(kappa_a + kappa_s))
    else
      cdt = clight * dt * 0.5 ! Half time step
    endif

    r2cdt2 = r**2 + cdt**2
    do j = -nmu+1, nmu-1
      r_new = sqrt(r2cdt2 + 2.0*mu(j)*r*cdt)
      mu_new = (mu(j)*r + cdt) / r_new

      ! Preserve direction, but rescale other directions
      rescale = sqrt( (1.0 - mu_new**2) / (u_com(2,1,j)**2 + u_com(3,1,j)**2) )

      ! Always transforming u_com, so all k are the same
      u_mu(0,:,j) = u_com(0,:,j)
      u_mu(1,:,j) = mu_new
      u_mu(2,:,j) = u_com(2,:,j) * rescale
      u_mu(3,:,j) = u_com(3,:,j) * rescale
    enddo

    u_mu(:,:,-nmu) = u_com(:,:,-nmu)
    u_mu(:,:, nmu) = u_com(:,:, nmu)

  end subroutine evolve_mu

  subroutine isotropic_core(f)
    use rice_config,    only: neps, nmu, npsi, nflav
    use rice_grid,      only: domega_scaled

    real, intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: fmean

    integer :: i, j, k, l

    do i = 1, neps
      do l = 1, nflav

        fmean = 0.0
        do j = -nmu, nmu
          do k = 1, npsi
            fmean = fmean + f(l,i,k,j) * domega_scaled(k,j)
          enddo
        enddo
        f(l,i,:,:) = fmean

      enddo
    enddo


  end subroutine isotropic_core

  subroutine isotropic_core_3d(f)
    use rice_config,    only: neps, nmu, npsi, nflav, nphi, ntheta, mpi_identical
    use rice_grid,      only: domega_scaled, volume, a_s, a_e, b_s, b_e, c_s, c_e
    use rice_mpi,       only: id, send_mpi, recv_mpi, bcast_mpi, sum_mpi

    real, intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)

    real    :: fsum, send, recv
    real    :: ivtot
    integer :: b, c, i, j, k, l
    integer :: tag
    logical :: has_cell

    ivtot = 1.0 / sum(volume(1:nphi,1:ntheta,1))

    do j = -nmu, nmu
      do k = 1, npsi

        do i = 1, neps
          do l = 1, nflav

            fsum = 0.0

            if (mpi_identical) then

              do b = 1, ntheta
                do c = 1, nphi

                  has_cell = .false.
                  if (a_s == 1) then ! Only send from tasks that have the center cell
                    if (b >= b_s .and. b <= b_e) then
                      if (c >= c_s .and. c <= c_e) then
                        has_cell = .true.
                      endif
                    endif
                  endif

                  if (has_cell) then
                    send = f(l,i,k,j,c,b,1) * domega_scaled(k,j) * volume(c,b,1)
                  endif

                  ! Use tag to ensure that master is receiving in the right order
                  tag = 10000*b+c

                  if (id == 0) then
                    if (has_cell) then
                      recv = send
                    else
                      call recv_mpi(recv, tag)
                    endif
                    fsum = fsum + recv
                  else
                    if (has_cell) then
                      call send_mpi(send, 0, tag)
                    endif
                  endif

                enddo
              enddo

              if (id == 0) then
                fsum = fsum * ivtot
              endif

              call bcast_mpi(fsum)

            else ! not mpi_identical

              if (a_s == 1) then
                do b = b_s, b_e
                  do c = c_s, c_e

                    ! distribution function * phase space volume/4pi * volume
                    fsum = fsum + f(l,i,k,j,c,b,1) * domega_scaled(k,j) * volume(c,b,1)

                  enddo
                enddo
                fsum = fsum * ivtot
              endif

              ! Sum across MPI tasks
              call sum_mpi(fsum)

            endif

            f(l,i,k,j,c_s:c_e,b_s:b_e,1) = fsum

          enddo
        enddo

      enddo
    enddo

  end subroutine isotropic_core_3d

  subroutine forward_core(f)
    use rice_config,    only: neps, nmu, npsi, nflav
    use rice_grid,      only: domega

    real,    intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: fmean

    integer :: i, j, k, l

    do i = 1, neps
      do l = 1, nflav

        fmean = 0.0

        do j = -nmu, 0
          do k = 1, npsi
            fmean = fmean + f(l,i,k,j) * domega(k,j)
          enddo
        enddo

        fmean = fmean / sum(domega(:,1:nmu))

        do j = -nmu, 0
          f(l,i,:,j) = 0.0
        enddo
        do j = 1, nmu
          f(l,i,:,j) = f(l,i,:,j) + fmean
        enddo

      enddo
    enddo

  end subroutine forward_core

  subroutine average_core(f)
    use rice_config,    only: neps, nmu, npsi, nflav

    real,    intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: fmean

    integer :: i, j, k, l

    do i = 1, neps
      do l = 1, nflav

        do j = 1, nmu
          do k = 1, npsi
            fmean = 0.5 * (f(l,i,k,j) + f(l,i,k,-j))
            f(l,i,k, j) = fmean
            f(l,i,k,-j) = fmean
          enddo
        enddo

      enddo
    enddo

  end subroutine average_core

  subroutine scatter_core(f)
    use rice_config,    only: neps, nmu, npsi, nflav
    use rice_grid,      only: domega

    real,    intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: fmean
    real :: denom
    integer :: i, j, k, l

    do i = 1, neps
      do l = 1, nflav

        fmean = 0.0
        denom = 0.0
        do j = -nmu, -nmu+1
          do k = 1, npsi
            fmean = fmean + f(l,i,k,-nmu) * domega(k,-nmu)
            denom = denom + domega(k,-nmu)
          enddo
        enddo

        fmean = fmean / denom

        do j = -nmu, -nmu+1
          f(l,i,:,j) = fmean
        enddo

      enddo
    enddo

  end subroutine scatter_core

  subroutine average_pole(f)
    use rice_config,    only: neps, nmu, npsi, nflav

    real,    intent(inout) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: fmean
    real :: inpsi

    integer :: i, j, l

    inpsi = 1.0 / real(npsi)

    do i = 1, neps
      do l = 1, nflav
        do j = -nmu, nmu

          fmean = sum(f(l,i,:,j)) * inpsi
          f(l,i,:,j) = fmean

        enddo
      enddo
    enddo

  end subroutine average_pole

end module rice_treatments
