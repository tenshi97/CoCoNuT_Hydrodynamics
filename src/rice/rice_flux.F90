module rice_flux
 implicit none

contains

  function flux(f, u, dtslope) result(adv_flux)
    use rice_config, only: neps, nmu, npsi, nflav, clight

    real,           intent(in)  :: f      (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,           intent(in)  :: u      (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, optional, intent(in)  :: dtslope(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: adv_flux(1:nflav,1:neps,1:npsi,-nmu:nmu)

    if (present(dtslope)) then
      adv_flux = clight * u * (f - u*dtslope)
    else
      adv_flux = clight * u * f
    endif

  end function flux

  function flux_mono(f, u) result(adv_flux)
    use rice_config, only: neps, nmu, npsi, nflav, clight

    real,           intent(in)  :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,           intent(in)  :: u(               1:npsi,-nmu:nmu)

    real :: adv_flux(1:nflav,1:neps,1:npsi,-nmu:nmu)

    integer :: j, k

    do j = -nmu, nmu
      do k = 1, npsi
        adv_flux(:,:,k,j) = clight * u(k,j) * f(:,:,k,j)
      enddo
    enddo

  end function flux_mono

  subroutine normal_velocity(u_norm, u, sin_row_theta, cos_row_theta, beta, direction)
    use rice_config, only: nmu, npsi, cartoon_grid
    real,    intent(out) :: u_norm(1:npsi,-nmu:nmu)
    real,    intent(in)  :: u(0:3,1:npsi,-nmu:nmu)
    real,    intent(in)  :: sin_row_theta
    real,    intent(in)  :: cos_row_theta
    real,    intent(in)  :: beta
    integer, intent(in)  :: direction

    if (direction == 1) then
      u_norm = u(1,:,:) ! no rotation for radial interfaces
    elseif (direction == 2) then
      u_norm(:,:) = -u(1,:,:)*sin(beta) + u(2,:,:)*cos(beta)
    elseif (direction == 3) then
      if (cartoon_grid) then
        ! theta = pi/2
        u_norm(:,:) = -u(1,:,:)*sin(beta) + u(3,:,:)*cos(beta)
      else
        u_norm(:,:) = -u(1,:,:)*sin(beta)*sin_row_theta - u(2,:,:)*sin(beta)*cos_row_theta + u(3,:,:)*cos(beta)
      endif
    endif

  end subroutine normal_velocity

end module rice_flux
