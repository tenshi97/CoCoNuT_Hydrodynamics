module rice_transform
  implicit none

  real, parameter :: hs2 = 0.5*sqrt(2.0)

contains

  function transform_frame(f, u_dst, dx, sin_theta, cos_theta, direction, no_rotate) result(f_rot)
    use rice_constants,       only: pi
    use rice_config,          only: neps, nmu, npsi, nflav, transform_rotate
    use rice_grid,            only: domega, domega_inv, epsvol, epsvol_inv, cpsi, spsi, eps, mu, psi
    use rice_utils,           only: magnitude
    use rice_transform_utils, only: rotate_vector_theta, rotate_vector_phi, &
                                    get_nearest_index, get_nearest_index_reverse

    real,    intent(in) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,    intent(in) :: u_dst(0:3,1:npsi,-nmu:nmu) ! Boosted velocities
    real,    intent(in) :: dx                         ! coordinate grid rotation (geometric)
    real,    intent(in) :: sin_theta                  ! sin(theta) (used only for direction==3)
    real,    intent(in) :: cos_theta                  ! cos(theta) (used only for direction==3)
    integer, intent(in) :: direction                  ! 1=r, 2=theta, 3=phi
    logical, optional, intent(in) :: no_rotate
    real :: f_rot(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: beta ! rotation angle

    real :: mu_new, psi_new, eps_new
    real :: u_new(0:3)

    real :: sinbeta, cosbeta
    real :: cospsii, cospsia, cospsib
    real :: sinpsii, sinpsia, sinpsib

    real :: umag, umag_lat, umag_lat2

    real :: idenom, itotal

    real :: f_scr(1:nflav), scr3, scr4
    real :: wpsia, wpsib, wmua, wmub, wepsa, wepsb
    real :: om_inv_aa, om_inv_ab, om_inv_ba, om_inv_bb
    real :: om

    real :: scr_a, scr_b
    real :: scr_aa, scr_ab, scr_ba, scr_bb

    integer :: i, j, k, l
    integer :: ia, ib, ja, jb, ka, kb

    logical :: rotate

    f_rot = 0.0

    if (present(no_rotate)) then
      rotate = .false.
    elseif (.not. transform_rotate) then
      rotate = .false.
    elseif (direction == 1) then
      rotate = .false.
    else
      rotate = .true.
    endif

    ! Only rotate if not in radial direction
    if (rotate) then
      beta = - dx
      sinbeta = sin(beta)
      cosbeta = cos(beta)
    endif

    do j = -nmu, nmu
      do k = 1, npsi

        ! Perform geometric rotation on the destination angles
        if (.not. rotate) then
          ! No rotation for r
          u_new = u_dst(:,k,j)
        elseif (direction == 2) then
          call rotate_vector_theta(u_dst(:,k,j), u_new, sinbeta, cosbeta)
        elseif (direction == 3) then
          call rotate_vector_phi(u_dst(:,k,j), u_new, sinbeta, cosbeta, sin_theta, cos_theta)
        endif

        ! Magnitude of spatial component of u, and lateral component
        umag_lat2 = u_new(2)**2 + u_new(3)**2
        umag = sqrt(u_new(1)**2 + umag_lat2)
        umag_lat = sqrt(umag_lat2)

        mu_new  = u_new(1) / umag

        if (umag_lat > 0.0) then ! abs(mu) < 1.0
          cospsii = u_new(2) / umag_lat
          sinpsii = u_new(3) / umag_lat
        else ! mu = -1 or 1, psi does not matter, set it based on k
          cospsii = cpsi(k)
          sinpsii = spsi(k)
        endif

        ! Find nearest destination indices for mu and psi
        ! Search directly for mu
        call get_nearest_index(mu_new, mu, 2*nmu+1, j+nmu+1, 2.0, ja, jb)
        ja = ja - nmu - 1
        jb = jb - nmu - 1

        ! Search directly for psi
        psi_new = atan2(sinpsii, cospsii)
        call get_nearest_index(psi_new, psi, npsi, k, 2.0*pi, ka, kb)

        ! To bypass atan2, search using cos(psi) and sin(psi)
        ! Decide based on which one is monotonic
        ! if (sinpsii > hs2) then
        !   call get_nearest_index_reverse(cospsii, cpsi, npsi, k, 0.0, ka, kb)
        ! elseif (sinpsii < -hs2) then
        !   call get_nearest_index(cospsii, cpsi, npsi, k, 0.0, ka, kb)
        ! elseif (cospsii > hs2) then
        !   call get_nearest_index(sinpsii, spsi, npsi, k, 0.0, ka, kb)
        ! else ! (cospsii < -0.5*sqrt(2.0))
        !   call get_nearest_index_reverse(sinpsii, spsi, npsi, k, 0.0, ka, kb)
        ! endif

        ! Determine mu weights
        wmub = (mu_new - mu(ja)) / (mu(jb) - mu(ja))
        wmua = 1.0 - wmub

        ! Determine psi weights
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

        om = domega(k,j)
        om_inv_aa = domega_inv(ka, ja)
        om_inv_ab = domega_inv(ka, jb)
        om_inv_ba = domega_inv(kb, ja)
        om_inv_bb = domega_inv(kb, jb)

        scr_aa = wmua * wpsia * om_inv_aa
        scr_ab = wmub * wpsia * om_inv_ab
        scr_ba = wmua * wpsib * om_inv_ba
        scr_bb = wmub * wpsib * om_inv_bb

        ! Angles are the same for all energy bins
        ! Now, get energy weighting
        do i = 1, neps
          if (neps > 1) then
            ! Shifted energy bin
            eps_new = eps(i) * u_dst(0,k,j)
            call get_nearest_index(eps_new, eps, neps, i, 0.0, ia, ib)

            ! Extrapolate - energy conservative
            ! wepsb = (eps_new - eps(ia)) / (eps(ib) - eps(ia))
            ! wepsa = 1.0 - wepsb

            ! Don't shift beyond grid boundaries - number conservative
            if (eps_new < eps(ia)) then
              wepsa = 1.0
              wepsb = 0.0
            elseif (eps_new > eps(ib)) then
              wepsa = 0.0
              wepsb = 1.0
            else
              wepsb = (eps_new - eps(ia)) / (eps(ib) - eps(ia))
              wepsa = 1.0 - wepsb
            endif

          else ! If only one energy bin
            ia = 1
            ib = 1
            wepsa = 1.0
            wepsb = 0.0
          endif

          scr_a = wepsa * epsvol_inv(ia)
          scr_b = wepsb * epsvol_inv(ib)

          f_scr = f(1:nflav,i,k,j) * om * epsvol(i)

          ! Loop unrolled for optimisation
          do l = 1, nflav
             scr3 = scr_a * f_scr (l)
             f_rot(l,ia,ka,ja) = f_rot(l,ia,ka,ja) + &
                  scr_aa * scr3
             f_rot(l,ia,ka,jb) = f_rot(l,ia,ka,jb) + &
                  scr_ab * scr3
             f_rot(l,ia,kb,ja) = f_rot(l,ia,kb,ja) + &
                  scr_ba * scr3
             f_rot(l,ia,kb,jb) = f_rot(l,ia,kb,jb) + &
                  scr_bb * scr3
          enddo

          do l = 1, nflav
             scr4 = scr_b * f_scr (l)
             f_rot(l,ib,ka,ja) = f_rot(l,ib,ka,ja) + &
                  scr_aa * scr4
             f_rot(l,ib,ka,jb) = f_rot(l,ib,ka,jb) + &
                  scr_ab * scr4
             f_rot(l,ib,kb,ja) = f_rot(l,ib,kb,ja) + &
                  scr_ba * scr4
             f_rot(l,ib,kb,jb) = f_rot(l,ib,kb,jb) + &
                  scr_bb * scr4
          enddo

        enddo
      enddo
    enddo

  end function transform_frame

  function transform_frame_cartoon(f, b, ia) result(f_rot)
    use rice_config, only: nmu, npsi, neps, nflav
    use rice_grid,   only: cartoon_weight, cartoon_index

    real,    intent(in) :: f(1:nflav,1:neps,1:npsi,-nmu:nmu)
    integer, intent(in) :: b
    integer, intent(in) :: ia

    real :: wpsia, wpsib
    real :: f_rot(1:nflav,1:neps,1:npsi,-nmu:nmu)

    integer :: k, ka, kb

    f_rot = 0.0

    do k = 1, npsi

      ! Use precomputed indices
      ka = cartoon_index(1,k,ia,b)
      kb = cartoon_index(2,k,ia,b)

      ! Use precomputed values
      wpsia = cartoon_weight(k,ia,b)
      wpsib = 1.0 - wpsia

      ! Omega not needed, because psi bins are the same size
      f_rot(1:nflav,1:neps,ka,-nmu:nmu) = f_rot(1:nflav,1:neps,ka,-nmu:nmu) + &
                                  wpsia * f    (1:nflav,1:neps,k, -nmu:nmu)

      f_rot(1:nflav,1:neps,kb,-nmu:nmu) = f_rot(1:nflav,1:neps,kb,-nmu:nmu) + &
                                  wpsib * f    (1:nflav,1:neps,k, -nmu:nmu)

    enddo

  end function transform_frame_cartoon

end module rice_transform
