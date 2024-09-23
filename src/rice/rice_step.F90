module rice_step
  implicit none

contains

  subroutine step(dt)
    use rice_grid,        only: f, f_eq, kappa_a, kappa_s, &
                                r, theta, phi, r_if, theta_if, phi_if, &
                                volume, area_r, area_theta, area_phi, &
                                alpha, beta, phiconf, &
                                vfluid, jump, u_lab, &
                                a_s, b_s, c_s, a_e, b_e, c_e, na, nb, nc
    use rice_config,      only: neps, nmu, npsi, nflav, &
                                nr_nolateral, flux_r, flux_theta, flux_phi
    use rice_phase,       only: phase_step
    use rice_boundaries,  only: fill_boundary_ghosts
    use rice_mpi,         only: exchange

    real, intent(in)    :: dt

    real :: dfdt(1:nflav,1:neps,1:npsi,-nmu:nmu,c_s:c_e,b_s:b_e,a_s:a_e)

    integer :: a, b, c

    ! Perform MPI exchange first, so that information is available
    ! to fill hydro ghosts
    call exchange
    call fill_boundary_ghosts

    ! Calculate all u_lab (which are re-used multiple times)
    !$omp parallel do
    do a = a_s-2, a_e+2
      do b = b_s-2, b_e+2
        do c = c_s-2, c_e+2

          ! Don't compute corners because they don't contain valid values
          if ((a < a_s .or. a > a_e) .and. (b < b_s .or. b > b_e)) cycle
          if ((b < b_s .or. b > b_e) .and. (c < c_s .or. c > c_e)) cycle
          if ((a < a_s .or. a > a_e) .and. (c < c_s .or. c > c_e)) cycle

          call compute_ulab(u_lab(0:3,1:npsi,-nmu:nmu,c,b,a), &
                            vfluid(1:3,c,b,a), &
                            alpha(c,b,a), &
                            beta(1:3,c,b,a), &
                            phiconf(c,b,a))
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! Set flux to zero
    dfdt = 0.0

    ! Solve for advection in the r direction
    if (flux_r) then
      do b = b_s, b_e
        do c = c_s, c_e
          call solve_advection_stripe(dfdt    (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a_s  :a_e  ), &
                                      f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a_s-2:a_e+2), &
                                      f_eq    (1:nflav,1:neps,                c,b,a_s-2:a_e+2), &
                                      kappa_a (1:nflav,1:neps,                c,b,a_s-2:a_e+2), &
                                      kappa_s (1:nflav,1:neps,                c,b,a_s-2:a_e+2), &
                                      vfluid  (1:3,                           c,b,a_s-2:a_e+2), &
                                      alpha   (                               c,b,a_s-2:a_e+2), &
                                      beta    (1:3,                           c,b,a_s-2:a_e+2), &
                                      phiconf (                               c,b,a_s-2:a_e+2), &
                                      u_lab   (0:3,           1:npsi,-nmu:nmu,c,b,a_s-2:a_e+2), &
                                      na, a_s, a_e, dt, &
                                      r     (    a_s-2:a_e+2), &
                                      r_if  (    a_s-2:a_e+1), &
                                      volume(c,b,a_s-1:a_e+1), &
                                      area_r(c,b,a_s-1:a_e  ), &
                                      jump  (c,b,a_s-1:a_e  ), &
                                      -1.0, -1.0, -1.0, 1)
        enddo
      enddo
    endif

    ! Solve for advection in the theta direction
    if (flux_theta) then
      !$omp parallel do
      do a = max(nr_nolateral+1, a_s), a_e
        do c = c_s, c_e
          call solve_advection_stripe(dfdt    (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b_s  :b_e  ,a), &
                                      f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b_s-2:b_e+2,a), &
                                      f_eq    (1:nflav,1:neps,                c,b_s-2:b_e+2,a), &
                                      kappa_a (1:nflav,1:neps,                c,b_s-2:b_e+2,a), &
                                      kappa_s (1:nflav,1:neps,                c,b_s-2:b_e+2,a), &
                                      vfluid  (1:3,                           c,b_s-2:b_e+2,a), &
                                      alpha   (                               c,b_s-2:b_e+2,a), &
                                      beta    (1:3,                           c,b_s-2:b_e+2,a), &
                                      phiconf (                               c,b_s-2:b_e+2,a), &
                                      u_lab   (0:3,           1:npsi,-nmu:nmu,c,b_s-2:b_e+2,a), &
                                      nb, b_s, b_e, dt, &
                                      theta     (  b_s-2:b_e+2  ), &
                                      theta_if  (  b_s-2:b_e+1  ), &
                                      volume    (c,b_s-1:b_e+1,a), &
                                      area_theta(c,b_s-1:b_e,  a), &
                                      jump      (c,b_s-1:b_e,  a), &
                                      -1.0, r(a), r_if(a-1), 2)
        enddo
      enddo
      !$omp end parallel do
    endif

    ! Solve for advection in the phi direction
    if (flux_phi) then
      !$omp parallel do
      do a = max(nr_nolateral+1, a_s), a_e
        do b = b_s, b_e
          call solve_advection_stripe(dfdt    (1:nflav,1:neps,1:npsi,-nmu:nmu,c_s  :c_e  ,b,a), &
                                      f       (1:nflav,1:neps,1:npsi,-nmu:nmu,c_s-2:c_e+2,b,a), &
                                      f_eq    (1:nflav,1:neps,                c_s-2:c_e+2,b,a), &
                                      kappa_a (1:nflav,1:neps,                c_s-2:c_e+2,b,a), &
                                      kappa_s (1:nflav,1:neps,                c_s-2:c_e+2,b,a), &
                                      vfluid  (1:3,                           c_s-2:c_e+2,b,a), &
                                      alpha   (                               c_s-2:c_e+2,b,a), &
                                      beta    (1:3,                           c_s-2:c_e+2,b,a), &
                                      phiconf (                               c_s-2:c_e+2,b,a), &
                                      u_lab   (0:3,           1:npsi,-nmu:nmu,c_s-2:c_e+2,b,a), &
                                      nc, c_s, c_e, dt, &
                                      phi     (c_s-2:c_e+2    ), &
                                      phi_if  (c_s-2:c_e+1    ), &
                                      volume  (c_s-1:c_e+1,b,a), &
                                      area_phi(c_s-1:c_e,  b,a), &
                                      jump    (c_s-1:c_e,  b,a), &
                                      theta(b), r(a), r_if(a-1), 3)
        enddo
      enddo
      !$omp end parallel do
    endif

    ! With all the fluxes added up, solve for the full step
    !$omp parallel do
    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e
          call phase_step(f      (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a), &
                          dfdt   (1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a), &
                          f_eq   (1:nflav,1:neps,                c,b,a), &
                          kappa_a(1:nflav,1:neps,                c,b,a), &
                          kappa_s(1:nflav,1:neps,                c,b,a), &
                          dt)
        enddo
      enddo
    enddo
    !$omp end parallel do

    call core_treatments

  end subroutine step

  subroutine solve_advection_stripe(dfdt_stripe, f, f_eq, kappa_a, kappa_s, &
                                    vfluid, alpha, beta, phiconf, u_lab, &
                                    n, i_s, i_e, dt, x, x_if, &
                                    volume, area, &
                                    jump, row_theta, row_r, row_rif, direction)
    use rice_config,       only: neps, nmu, npsi, nflav, nr_diffusive, lw_rotate, cartoon_grid
    use rice_laxwendroff,  only: lw_state, interpolation_weight
    use rice_upwind,       only: upwind_state
    use rice_reconstruct,  only: reconstruct_linear
    use rice_transform,    only: transform_frame
    use rice_boost,        only: boost, boost_split
    use rice_grid,         only: u_com, f_null
    use rice_sr,           only: lambda_transform
    use rice_gr,           only: m_transform, m_inv_transform, metric, metric_inv, &
                                 ricci_rotation, gr_flux_correction, gr_flux_u0_divide
    use rice_flux,         only: flux, normal_velocity
    use rice_treatments,   only: evolve_mu

    integer, intent(in)     :: n
    integer, intent(in)     :: i_s
    integer, intent(in)     :: i_e
    real,    intent(inout)  :: dfdt_stripe(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s  :i_e  ) ! calculated flux
    real,    intent(in)     :: f          (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-2:i_e+2)
    real,    intent(in)     :: f_eq       (1:nflav,1:neps,                i_s-2:i_e+2)
    real,    intent(in)     :: kappa_a    (1:nflav,1:neps,                i_s-2:i_e+2)
    real,    intent(in)     :: kappa_s    (1:nflav,1:neps,                i_s-2:i_e+2)
    real,    intent(in)     :: vfluid     (1:3,                           i_s-2:i_e+2)
    real,    intent(in)     :: alpha      (                               i_s-2:i_e+2)
    real,    intent(in)     :: beta       (1:3,                           i_s-2:i_e+2)
    real,    intent(in)     :: phiconf    (                               i_s-2:i_e+2)
    real,    intent(in)     :: u_lab      (0:3,           1:npsi,-nmu:nmu,i_s-2:i_e+2)
    real,    intent(in)     :: dt
    real,    intent(in)     :: x     (i_s-2:i_e+2) ! cell centers
    real,    intent(in)     :: x_if  (i_s-2:i_e+1) ! cell interfaces
    real,    intent(in)     :: volume(i_s-1:i_e+1)
    real,    intent(in)     :: area  (i_s-1:i_e)
    logical, intent(in)     :: jump  (i_s-1:i_e)
    real,    intent(in)     :: row_theta ! theta (used if solving along phi)
    real,    intent(in)     :: row_r     ! r (used if solving along lateral directions)
    real,    intent(in)     :: row_rif   ! r_if (used for flux switching)
    integer, intent(in)     :: direction ! 1=r, 2=theta, 3=phi

    ! Fluid velocity at interface
    real :: vfluid_if(3)

    ! 3+1 CFC terms at interfaces
    real :: alpha_if
    real :: beta_if(3)
    real :: phiconf_if

    ! Lorentz boost matrices
    real :: lambda_l     (0:3,0:3)
    real :: lambda_r     (0:3,0:3)
    real :: lambda_if    (0:3,0:3)
    real :: lambda_inv_l (0:3,0:3)
    real :: lambda_inv_r (0:3,0:3)
    real :: lambda_inv_if(0:3,0:3)

    ! Eulerian -> coordinate boost matrices
    real :: m_l     (0:3,0:3)
    real :: m_r     (0:3,0:3)
    real :: m_if    (0:3,0:3)
    real :: m_inv_l (0:3,0:3)
    real :: m_inv_r (0:3,0:3)
    real :: m_inv_if(0:3,0:3)

    ! Full boost matrices
    real :: mfull_inv_l (0:3,0:3)
    real :: mfull_inv_r (0:3,0:3)
    real :: mfull_inv_if(0:3,0:3)

    ! Metric
    real :: g_l     (0:3,0:3)
    real :: g_r     (0:3,0:3)
    real :: g_if    (0:3,0:3)
    real :: g_inv_l (0:3,0:3)
    real :: g_inv_r (0:3,0:3)
    real :: g_inv_if(0:3,0:3)

    ! Fluxes
    real :: flux_l        (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_r        (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_upw_l_src(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_upw_r_src(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_upw_l_dst(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_upw_r_dst(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_lw       (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_lw_l     (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: flux_lw_r     (1:nflav,1:neps,1:npsi,-nmu:nmu)

    ! Intermediate distributions
    real :: f_upw_l  (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: f_upw_r  (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: f_lw     (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: f_l_if   (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: f_r_if   (1:nflav,1:neps,1:npsi,-nmu:nmu)

    ! Opacity switched velocities
    real :: u_if_lab_weighted     (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: u_l_lab_weighted      (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: u_r_lab_weighted      (1:nflav,1:neps,1:npsi,-nmu:nmu)

    ! Destination bins for phase space
    real :: u_l_r         (0:3,1:npsi,-nmu:nmu)
    real :: u_r_l         (0:3,1:npsi,-nmu:nmu)
    real :: u_if_l        (0:3,1:npsi,-nmu:nmu)
    real :: u_if_r        (0:3,1:npsi,-nmu:nmu)
    real :: u_l_if        (0:3,1:npsi,-nmu:nmu)
    real :: u_r_if        (0:3,1:npsi,-nmu:nmu)

    ! Intermediate advection velocities
    real :: u_adv_if_eul  (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_if_eul(0:3,1:npsi,-nmu:nmu)
    real :: u_adv_l_eul   (0:3,1:npsi,-nmu:nmu)
    real :: u_adv_r_eul   (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_l_eul (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_r_eul (0:3,1:npsi,-nmu:nmu)
    real :: u_adv_if_lab  (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_if_lab(0:3,1:npsi,-nmu:nmu)
    real :: u_if_l_lab    (0:3,1:npsi,-nmu:nmu)
    real :: u_if_r_lab    (0:3,1:npsi,-nmu:nmu)
    real :: u_l_r_lab     (0:3,1:npsi,-nmu:nmu)
    real :: u_r_l_lab     (0:3,1:npsi,-nmu:nmu)
    real :: u_adv_l_lab   (0:3,1:npsi,-nmu:nmu)
    real :: u_adv_r_lab   (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_l_lab (0:3,1:npsi,-nmu:nmu)
    real :: u_trans_r_lab (0:3,1:npsi,-nmu:nmu)
    real :: u_l_lab       (0:3,1:npsi,-nmu:nmu)
    real :: u_r_lab       (0:3,1:npsi,-nmu:nmu)
    real :: u_if_lab      (0:3,1:npsi,-nmu:nmu)

    ! Single component velocities
    real :: u_adv_l_lab_norm  (1:npsi,-nmu:nmu)
    real :: u_adv_r_lab_norm  (1:npsi,-nmu:nmu)
    real :: u_trans_l_lab_norm(1:npsi,-nmu:nmu)
    real :: u_trans_r_lab_norm(1:npsi,-nmu:nmu)

    ! Linear reconstruction
    real :: slope_l  (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real :: slope_r  (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real :: f_recon_l(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real :: f_recon_r(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)

    ! dfdt
    real :: dfdt_stripe_l(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s:i_e)
    real :: dfdt_stripe_r(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s:i_e)

    ! Precomputed trig values
    real :: sin_row_theta, cos_row_theta

    ! Switch weights
    real :: w_l, w_r

    real :: r_l, r_r
    real :: dx
    integer :: i

    logical :: omp_in_parallel

    ! Precompute sine and cos
    if (direction == 3) then
      sin_row_theta = sin(row_theta)
      cos_row_theta = cos(row_theta)
    endif

    call reconstruct_linear(n, i_s, i_e, f, x, x_if, slope_l, slope_r, f_recon_l, f_recon_r, vfluid, &
                            alpha, beta, phiconf, u_lab, sin_row_theta, cos_row_theta, direction)

    ! Over interfaces
    !$omp parallel if (i_e - i_s + 1 > 1 .and. .not. omp_in_parallel()) &
    !$omp default(none) &
    !$omp private(i, f_upw_l, f_upw_r, f_lw) &
    !$omp private(u_adv_l_eul, u_adv_r_eul, u_trans_l_eul, u_trans_r_eul) &
    !$omp private(u_adv_if_eul, u_trans_if_eul) &
    !$omp private(u_adv_l_lab, u_adv_r_lab, u_trans_l_lab, u_trans_r_lab) &
    !$omp private(u_adv_if_lab, u_trans_if_lab, u_l_lab_weighted, u_r_lab_weighted) &
    !$omp private(u_if_lab_weighted) &
    !$omp private(u_l_lab, u_r_lab, u_if_lab) &
    !$omp private(u_adv_l_lab_norm, u_trans_l_lab_norm, u_adv_r_lab_norm, u_trans_r_lab_norm) &
    !$omp private(flux_l, flux_r, flux_upw_l_src, flux_upw_r_src, flux_upw_l_dst, flux_upw_r_dst) &
    !$omp private(flux_lw, flux_lw_l, flux_lw_r) &
    !$omp private(u_r_l_lab, u_l_r_lab, u_if_l_lab, u_if_r_lab) &
    !$omp private(u_l_r, u_r_l, u_if_l, u_if_r) &
    !$omp private(u_l_if, u_r_if, f_l_if, f_r_if) &
    !$omp private(vfluid_if, w_l, w_r) &
    !$omp private(r_l, r_r, dx) &
    !$omp private(alpha_if, beta_if, phiconf_if) &
    !$omp private(m_l, m_r, m_if, m_inv_l, m_inv_r, m_inv_if) &
    !$omp private(lambda_l, lambda_r, lambda_if, lambda_inv_l, lambda_inv_r, lambda_inv_if) &
    !$omp private(mfull_inv_l, mfull_inv_r, mfull_inv_if) &
    !$omp private(g_l, g_r, g_if, g_inv_l, g_inv_r, g_inv_if) &
    !$omp shared(u_com) &
    !$omp shared(nr_diffusive) &
    !$omp shared(n, i_s, i_e, f, vfluid, x_if, x, direction, sin_row_theta, cos_row_theta) &
    !$omp shared(f_eq, kappa_a, kappa_s, f_null, dt, area, volume) &
    !$omp shared(slope_l, slope_r, f_recon_l, f_recon_r) &
    !$omp shared(dfdt_stripe_l, dfdt_stripe_r) &
    !$omp shared(jump, lw_rotate, row_r, row_rif, cartoon_grid) &
    !$omp shared(u_lab) &
    !$omp shared(alpha, beta, phiconf)
    !$omp do
    do i = i_s-1, i_e
      if (direction == 1) then
        r_l = x(i  )
        r_r = x(i+1)
      else
        r_l = row_r
        r_r = row_r
      endif

      ! Get interpolated velocity
      call interpolation_weight(x(i), x(i+1), x_if(i), w_l, w_r)
      vfluid_if = vfluid(:,i)*w_l + vfluid(:,i+1)*w_r

      ! Interface lapse, shift, and conformal factor
      alpha_if   = alpha(i)   * w_l + alpha(i+1)   * w_r
      beta_if(:) = beta(:,i)  * w_l + beta(:,i+1)  * w_r
      phiconf_if = phiconf(i) * w_l + phiconf(i+1) * w_r

      !---------------------------------------------------------------------------------------------------------------------------
      ! Boosts for Upwind
      !

      ! Metrics
      g_l      = metric    (alpha(i),   beta(:,i),   phiconf(i)  )
      g_r      = metric    (alpha(i+1), beta(:,i+1), phiconf(i+1))
      g_if     = metric    (alpha_if,   beta_if(:),  phiconf_if  )
      g_inv_l  = metric_inv(alpha(i),   beta(:,i),   phiconf(i)  )
      g_inv_r  = metric_inv(alpha(i+1), beta(:,i+1), phiconf(i+1))
      g_inv_if = metric_inv(alpha_if,   beta_if(:),  phiconf_if  )

      ! SR boost matrices
      lambda_l      = lambda_transform(vfluid(:,i  ) )
      lambda_r      = lambda_transform(vfluid(:,i+1) )
      lambda_if     = lambda_transform(vfluid_if     )
      lambda_inv_l  = lambda_transform(-vfluid(:,i  ))
      lambda_inv_r  = lambda_transform(-vfluid(:,i+1))
      lambda_inv_if = lambda_transform(-vfluid_if    )

      ! GR boost matrices
      m_l      = m_transform    (alpha(i),   beta(:,i),   phiconf(i)  )
      m_r      = m_transform    (alpha(i+1), beta(:,i+1), phiconf(i+1))
      m_if     = m_transform    (alpha_if,   beta_if,     phiconf_if  )
      m_inv_l  = m_inv_transform(alpha(i),   beta(:,i),   phiconf(i)  )
      m_inv_r  = m_inv_transform(alpha(i+1), beta(:,i+1), phiconf(i+1))
      m_inv_if = m_inv_transform(alpha_if,   beta_if,     phiconf_if  )

      ! Full boosts
      mfull_inv_l  = matmul(lambda_inv_l,  m_inv_l )
      mfull_inv_r  = matmul(lambda_inv_r,  m_inv_r )
      mfull_inv_if = matmul(lambda_inv_if, m_inv_if)

      ! Boost from source to lab frame to get advection velocity
      ! Left to lab
      call boost_split(lambda_l, u_com, u_adv_l_eul, u_trans_l_eul)
      call boost(m_l, u_adv_l_eul,   u_adv_l_lab)
      call boost(m_l, u_trans_l_eul, u_trans_l_lab)
      u_l_lab = u_lab(:,:,:,i) ! same as u_l_lab = u_adv_l_lab + u_trans_l_lab

      ! Right to lab
      call boost_split(lambda_r, u_com, u_adv_r_eul, u_trans_r_eul)
      call boost(m_r, u_adv_r_eul,   u_adv_r_lab)
      call boost(m_r, u_trans_r_eul, u_trans_r_lab)
      u_r_lab = u_lab(:,:,:,i+1) ! Same as: u_r_lab = u_adv_r_lab + u_trans_r_lab

      ! Get the velocity through the interface as a result of the rotation
      call normal_velocity(u_adv_l_lab_norm,   u_adv_l_lab,   sin_row_theta, cos_row_theta, x_if(i) - x(i),   direction)
      call normal_velocity(u_trans_l_lab_norm, u_trans_l_lab, sin_row_theta, cos_row_theta, x_if(i) - x(i),   direction)
      call normal_velocity(u_adv_r_lab_norm,   u_adv_r_lab,   sin_row_theta, cos_row_theta, x_if(i) - x(i+1), direction)
      call normal_velocity(u_trans_r_lab_norm, u_trans_r_lab, sin_row_theta, cos_row_theta, x_if(i) - x(i+1), direction)

      !---------------------------------------------------------------------------------------------------------------------------
      ! Boosts for Lax-Wendroff
      !

      ! Boost from the left and right frames to the interface frame for LW interpolation
      ! A straight-line interpolation is performed (i.e. energy boost only)
      ! Bypass kick and approximate by going straight to next cell
      call boost(mfull_inv_if, u_l_lab, u_l_if)
      call boost(mfull_inv_if, u_r_lab, u_r_if)

      if (lw_rotate) then
        f_l_if = transform_frame(f(:,:,:,:,i)  , u_l_if, x_if(i) - x(i  ), sin_row_theta, cos_row_theta, direction)
        f_r_if = transform_frame(f(:,:,:,:,i+1), u_r_if, x_if(i) - x(i+1), sin_row_theta, cos_row_theta, direction)
      else
        f_l_if = transform_frame(f(:,:,:,:,i  ), u_l_if, 0.0, 0.0, 0.0, direction, no_rotate=.true.)
        f_r_if = transform_frame(f(:,:,:,:,i+1), u_r_if, 0.0, 0.0, 0.0, direction, no_rotate=.true.)
      endif

      ! Boost from the interface frame to the lab frame to get advection velocity
      call boost_split(lambda_if, u_com, u_adv_if_eul, u_trans_if_eul)
      call boost(m_if, u_adv_if_eul,   u_adv_if_lab  )
      call boost(m_if, u_trans_if_eul, u_trans_if_lab)
      u_if_lab = u_adv_if_lab + u_trans_if_lab

      !---------------------------------------------------------------------------------------------------------------------------
      ! Switch between LW flux and upwind flux based on kappa using the advection velocity
      !

      if (direction == 1) then
        dx = x(i+1) - x(i)
      else
        dx = 4.0 * (row_r - row_rif)
      endif
      call advection_velocity_weight(u_l_lab_weighted, u_r_lab_weighted, u_if_lab_weighted, &
                                     u_adv_l_lab_norm, u_adv_r_lab_norm, &
                                     u_trans_l_lab_norm, u_trans_r_lab_norm, &
                                     u_trans_if_lab, &
                                     kappa_a(1:nflav,1:neps,i), kappa_a(1:nflav,1:neps,i+1), &
                                     kappa_s(1:nflav,1:neps,i), kappa_s(1:nflav,1:neps,i+1), &
                                     jump(i), dx, direction)

      !---------------------------------------------------------------------------------------------------------------------------
      ! Upwind
      !

      ! Get the upwind part of the distribution travelling from left to right
      call upwind_state(f_upw_l, &
                        f_recon_r(:,:,:,:,i), &
                        f_null, &
                        u_l_lab_weighted)

      ! Get the upwind part of the distribution travelling from right to left
      call upwind_state(f_upw_r, &
                        f_null, &
                        f_recon_l(:,:,:,:,i+1), &
                        u_r_lab_weighted)

      !---------------------------------------------------------------------------------------------------------------------------
      ! Lax-Wendroff
      !

      ! Get the Lax-Wendroff state at the interface
      call lw_state(f_lw, &
                    f_l_if, &
                    f_r_if, &
                    x_if(i), &
                    x(i), &
                    x(i+1), &
                    f_eq(1:nflav,1:neps,i  ), &
                    f_eq(1:nflav,1:neps,i+1), &
                    kappa_a(1:nflav,1:neps,i  ), &
                    kappa_a(1:nflav,1:neps,i+1), &
                    kappa_s(1:nflav,1:neps,i  ), &
                    kappa_s(1:nflav,1:neps,i+1), &
                    volume(i), volume(i+1), area(i), &
                    dt, &
                    (direction == 1 .and. i < nr_diffusive), &
                    direction)

      !---------------------------------------------------------------------------------------------------------------------------
      ! Combine fluxes
      !

      ! Upwind flux sourced from left (no GR correction for source side)
      flux_upw_l_src = flux(f_upw_l, u_l_lab_weighted, dt*slope_r(:,:,:,:,i))

      ! Upwind flux sourced from right (no GR correction for source side)
      flux_upw_r_src = flux(f_upw_r, u_r_lab_weighted, dt*slope_l(:,:,:,:,i+1))

      ! Upwind arriving at the left

      ! u_adv_r_lab + u_trans_r_lab
      ! calculated above is sourced from the right frame, boosted to lab frame
      ! Now, boost into left frame to get the destination bin
      call ricci_rotation(g_r, g_inv_r, g_inv_l, u_r_lab, u_r_l_lab, direction)
      call boost(mfull_inv_l, u_r_l_lab, u_r_l)
      flux_upw_l_dst = transform_frame(flux_upw_r_src, u_r_l, x(i) - x(i+1), sin_row_theta, cos_row_theta, direction)

      ! Upwind arriving at the right
      call ricci_rotation(g_l, g_inv_l, g_inv_r, u_l_lab, u_l_r_lab, direction)
      call boost(mfull_inv_r, u_l_r_lab, u_l_r)
      flux_upw_r_dst = transform_frame(flux_upw_l_src, u_l_r, x(i+1) - x(i), sin_row_theta, cos_row_theta, direction)

      flux_lw = flux(f_lw, u_if_lab_weighted)

      ! LW sourced from left to interface and sourced from interface to left
      call ricci_rotation(g_if, g_inv_if, g_inv_l, u_if_lab, u_if_l_lab, direction)
      call boost(mfull_inv_l, u_if_l_lab, u_if_l)
      if (lw_rotate) then
       flux_lw_l = transform_frame(flux_lw, u_if_l, x(i) - x_if(i), sin_row_theta, cos_row_theta, direction)
      else
       flux_lw_l = transform_frame(flux_lw, u_if_l, 0.0, 0.0, 0.0, direction)
      endif

       ! LW sourced from right to interface and sourced from interface to right
      call ricci_rotation(g_if, g_inv_if, g_inv_r, u_if_lab, u_if_r_lab, direction)
      call boost(mfull_inv_r, u_if_r_lab, u_if_r)
      if (lw_rotate) then
       flux_lw_r = transform_frame(flux_lw, u_if_r, x(i+1) - x_if(i), sin_row_theta, cos_row_theta, direction)
      else
       flux_lw_r = transform_frame(flux_lw, u_if_r, 0.0, 0.0, 0.0, direction)
      endif

      ! Volume and area corrections
      call gr_flux_correction(flux_upw_l_dst, phiconf(i),   phiconf(i+1), alpha(i),   alpha(i+1))
      call gr_flux_correction(flux_upw_r_dst, phiconf(i+1), phiconf(i),   alpha(i+1), alpha(i)  )
      call gr_flux_correction(flux_lw_l,      phiconf(i)  , phiconf_if,   alpha(i),   alpha_if  )
      call gr_flux_correction(flux_lw_r,      phiconf(i+1), phiconf_if,   alpha(i+1), alpha_if  )

      ! Add up all components of flux
      flux_l = flux_upw_l_src + flux_upw_l_dst + flux_lw_l
      flux_r = flux_upw_r_src + flux_upw_r_dst + flux_lw_r

      ! Divide by u0
      call gr_flux_u0_divide(flux_l, u_l_lab)
      call gr_flux_u0_divide(flux_r, u_r_lab)

      ! Add left and right fluxes
      if (i   >= i_s) dfdt_stripe_r(:,:,:,:,i)   = - flux_l * area(i) / volume(i)
      if (i+1 <= i_e) dfdt_stripe_l(:,:,:,:,i+1) =   flux_r * area(i) / volume(i+1)
    enddo
    !$omp enddo
    !$omp end parallel

    ! Rather than doing an omp reduction on dfdt_stripe,
    ! save the left and right stripes and combine afterwards
    ! because each cell is only accessed twice
    ! stripe_l is sourced from the left interface
    dfdt_stripe = dfdt_stripe + dfdt_stripe_l + dfdt_stripe_r

  end subroutine solve_advection_stripe

  subroutine advection_velocity_weight(u_l_lab, u_r_lab, u_if_lab, &
                                       u_adv_l_lab, u_adv_r_lab, u_trans_l_lab, u_trans_r_lab, u_trans_if_lab, &
                                       kappa_a_l, kappa_a_r, kappa_s_l, kappa_s_r, &
                                       jump, dx, direction)
    use rice_config, only: neps, nmu, npsi, nflav, lw_only

    real, intent(out)  :: u_l_lab (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(out)  :: u_r_lab (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(out)  :: u_if_lab(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real, intent(in)   :: u_adv_l_lab       (1:npsi,-nmu:nmu)
    real, intent(in)   :: u_adv_r_lab       (1:npsi,-nmu:nmu)
    real, intent(in)   :: u_trans_l_lab     (1:npsi,-nmu:nmu)
    real, intent(in)   :: u_trans_r_lab     (1:npsi,-nmu:nmu)
    real, intent(in)   :: u_trans_if_lab(0:3,1:npsi,-nmu:nmu)

    real, intent(in)   :: kappa_a_l(1:nflav,1:neps)
    real, intent(in)   :: kappa_a_r(1:nflav,1:neps)
    real, intent(in)   :: kappa_s_l(1:nflav,1:neps)
    real, intent(in)   :: kappa_s_r(1:nflav,1:neps)

    logical, intent(in) :: jump
    real,    intent(in) :: dx
    integer, intent(in) :: direction

    real :: tau
    real :: wsplit
    real :: kappa

    integer :: i, l

    do i = 1, neps
      do l = 1, nflav
        if (jump) then
          wsplit = 1.0
          kappa = 0.0
          tau = 0.0
        else
          ! kappa = min(kappa_a_l(l,i) + kappa_s_l(l,i), kappa_a_r(l,i) + kappa_s_r(l,i))
          kappa = 0.5 * (kappa_a_l(l,i) + kappa_s_l(l,i) + kappa_a_r(l,i) + kappa_s_r(l,i))
          ! Weight based on opacity
          tau = kappa * dx
          wsplit = exp(-tau)
        endif

        ! Always upwind
        ! wsplit = 1.0

        if (lw_only) then
          wsplit = 0.0
        endif

        ! Upwind velocity
        u_l_lab(l,i,:,:) = u_adv_l_lab + wsplit * u_trans_l_lab
        u_r_lab(l,i,:,:) = u_adv_r_lab + wsplit * u_trans_r_lab

        ! LW velocity
        u_if_lab(l,i,:,:) = (1.0 - wsplit) * u_trans_if_lab(direction,:,:)
      enddo
    enddo

  end subroutine advection_velocity_weight

  subroutine time_boost
    use rice_config,    only: nmu, npsi, neps, nflav
    use rice_grid,      only: f, u_com, a_s, a_e, b_s, b_e, c_s, c_e, &
                              vfluid_old, vfluid, alpha_old, alpha, &
                              beta_old, beta, phiconf_old, phiconf
    use rice_boost,     only: boost
    use rice_sr,        only: lambda_transform
    use rice_gr,        only: metric, metric_inv, &
                              m_transform, m_inv_transform, ricci_rotation
    use rice_transform, only: transform_frame

    real :: u_old_lab(0:3,1:npsi,-nmu:nmu)
    real :: u_new_lab(0:3,1:npsi,-nmu:nmu)
    real :: u_new    (0:3,1:npsi,-nmu:nmu)

    real :: f_new(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: g_old         (0:3,0:3)
    real :: g_inv_old     (0:3,0:3)
    real :: g_inv_new     (0:3,0:3)
    real :: m_old         (0:3,0:3)
    real :: m_inv_new     (0:3,0:3)
    real :: lambda_old    (0:3,0:3)
    real :: lambda_inv_new(0:3,0:3)

    integer :: a, b, c

    !$omp parallel &
    !$omp default(none) &
    !$omp private(a, b, c) &
    !$omp private(g_old, g_inv_old, g_inv_new) &
    !$omp private(m_old, m_inv_new) &
    !$omp private(lambda_old, lambda_inv_new) &
    !$omp private(u_old_lab, u_new_lab, u_new) &
    !$omp private(f_new) &
    !$omp shared(a_s, a_e, b_s, b_e, c_s, c_e) &
    !$omp shared(vfluid_old, vfluid) &
    !$omp shared(phiconf_old, phiconf, alpha_old, alpha, beta_old, beta) &
    !$omp shared(u_com) &
    !$omp shared(f)
    !$omp do

    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e

          ! Metrics
          g_old     = metric    (alpha_old(c,b,a), beta_old(:,c,b,a), phiconf_old(c,b,a))
          g_inv_old = metric_inv(alpha_old(c,b,a), beta_old(:,c,b,a), phiconf_old(c,b,a))
          g_inv_new = metric_inv(alpha    (c,b,a), beta    (:,c,b,a), phiconf    (c,b,a))

          ! SR boost matrices
          lambda_old     = lambda_transform( vfluid_old(:,c,b,a))
          lambda_inv_new = lambda_transform(-vfluid    (:,c,b,a))

          ! GR boost matrices
          m_old     = m_transform    (alpha_old(c,b,a), beta_old(:,c,b,a), phiconf_old(c,b,a))
          m_inv_new = m_inv_transform(alpha    (c,b,a), beta    (:,c,b,a), phiconf    (c,b,a))

          ! Boost to lab frame
          call boost(matmul(m_old, lambda_old), u_com, u_old_lab)

          ! Boost to new time step in lab frame
          call ricci_rotation(g_old, g_inv_old, g_inv_new, u_old_lab, u_new_lab, 0)

          ! Boost to comoving frame
          call boost(matmul(lambda_inv_new, m_inv_new), u_new_lab, u_new)

          ! Transform frame
          f_new = transform_frame(f(:,:,:,:,c,b,a), u_new, 0.0, 0.0, 0.0, 0, no_rotate=.true.)

          ! Account for change in volume
          f(:,:,:,:,c,b,a) = f_new * (phiconf_old(c,b,a) / phiconf(c,b,a)) ** 6
        enddo
      enddo
    enddo
    !$omp enddo
    !$omp end parallel

    vfluid_old  = vfluid
    alpha_old   = alpha
    beta_old    = beta
    phiconf_old = phiconf

  end subroutine time_boost

  subroutine core_treatments
    use rice_config,      only: ntheta, neps, nmu, npsi, nflav, &
                                nr_isotropic, forward_inner, scatter_inner, &
                                average_inner, average_poles, nr_evolve_mu, &
                                isotropic_3d
    use rice_grid,        only: f, r, area_theta, kappa_a, kappa_s, &
                                a_s, a_e, b_s, b_e, c_s, c_e
    use rice_treatments,  only: evolve_mu, isotropic_core, forward_core, &
                                average_core, scatter_core, average_pole, &
                                isotropic_core_3d
    use rice_transform,   only: transform_frame

    real :: u_mu(0:3,1:npsi,-nmu:nmu)
    real :: f_new(1:nflav,1:neps,1:npsi,-nmu:nmu)

    integer :: a, b, c

    if (isotropic_3d) then
      ! - Works across MPI
      ! - Needs to happen before local treatments
      call isotropic_core_3d(f)
    endif

    do a = a_s, min(a_e, nr_isotropic)
      do b = b_s, b_e
        do c = c_S, c_e
          ! Make central zones isotropic
          call isotropic_core(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a))
        enddo
      enddo
    enddo

    if (forward_inner) then
      if ((a_s <= 1) .and. (a_e >= 1)) then
        do b = b_s, b_e
          do c = c_s, c_e
            ! Make central zones mu=1
            call forward_core(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,1))
          enddo
        enddo
      endif
    endif

    if (scatter_inner) then
      if ((a_s <= 1) .and. (a_e >= 1)) then
        do b = b_s, b_e
          do c = c_s, c_e
            ! Redistribute mu=-1
            call scatter_core(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,1))
          enddo
        enddo
      endif
    endif

    if (average_inner) then
      if ((a_s <= 1) .and. (a_e >= 1)) then
        do b = b_s, b_e
          do c = c_s, c_e
            ! Average -mu and +mu
            call average_core(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,1))
          enddo
        enddo
      endif
    endif

    if (average_poles) then
      if (area_theta(1,0,1) <= 0.0) then
        if ((b_s <= 1) .and. (b_e >= 1)) then
          do a = a_s, a_e
            do c = c_s, c_e
              call average_pole(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,1,a))
            enddo
          enddo
        endif
      endif
      if (area_theta(1,ntheta,1) <= 0.0) then
        if ((b_s <= ntheta) .and. (b_e >= ntheta)) then
          do a = a_s, a_e
            do c = c_s, c_e
              call average_pole(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,ntheta,a))
            enddo
          enddo
        endif
      endif
    endif

    do a = a_s, min(a_e, nr_evolve_mu)
      do b = b_s, b_e
        do c = c_s, c_e
          ! Evolve mu within cell
          call evolve_mu(r(a), 0.0, u_mu, kappa_a(1:nflav,1:neps,c,b,a), kappa_s(1:nflav,1:neps,c,b,a))
          f_new = transform_frame(f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a), u_mu, 0.0, 0.0, 0.0, 0, no_rotate=.true.)
          f(1:nflav,1:neps,1:npsi,-nmu:nmu,c,b,a) = f_new
        enddo
      enddo
    enddo

  end subroutine core_treatments

  subroutine compute_ulab(u_lab, vfluid, alpha, beta, phiconf)
    use rice_config,     only: nmu, npsi
    use rice_boost,      only: boost
    use rice_sr,         only: lambda_transform
    use rice_gr,         only: m_transform
    use rice_grid,       only: u_com
    real, intent(out) :: u_lab(0:3,1:npsi,-nmu:nmu)
    real, intent(in)  :: vfluid(1:3)
    real, intent(in)  :: alpha
    real, intent(in)  :: beta(1:3)
    real, intent(in)  :: phiconf

    call boost(matmul(m_transform(alpha, beta, phiconf), lambda_transform(vfluid)), &
               u_com, u_lab)

  end subroutine compute_ulab

end module rice_step
