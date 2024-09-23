module rice_reconstruct
  implicit none

contains

  subroutine reconstruct_linear(n, i_s, i_e, f, x, x_if, slope_l, slope_r, f_recon_l, f_recon_r, vfluid, &
                                alpha, beta, phiconf, u_lab, sin_row_theta, cos_row_theta, direction)
    use rice_config,     only: neps, nmu, npsi, nflav, reconstruct_fr2, slope_limiter
    use rice_transform,  only: transform_frame
    use rice_boost,      only: boost
    use rice_sr,         only: lambda_transform
    use rice_gr,         only: ricci_rotation, m_transform, m_inv_transform, &
                               metric, metric_inv

    integer,  intent(in)  :: n
    integer,  intent(in)  :: i_s
    integer,  intent(in)  :: i_e
    real,     intent(in)  :: f        (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-2:i_e+2)
    real,     intent(in)  :: x_if                                    (i_s-2:i_e+1)
    real,     intent(in)  :: x                                       (i_s-2:i_e+2)
    real,     intent(out) :: slope_l  (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real,     intent(out) :: slope_r  (1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real,     intent(out) :: f_recon_l(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real,     intent(out) :: f_recon_r(1:nflav,1:neps,1:npsi,-nmu:nmu,i_s-1:i_e+1)
    real,     intent(in)  :: vfluid                                (3,i_s-2:i_e+2)
    real,     intent(in)  :: alpha                                 (  i_s-2:i_e+2)
    real,     intent(in)  :: beta                                  (3,i_s-2:i_e+2)
    real,     intent(in)  :: phiconf                               (  i_s-2:i_e+2)
    real,     intent(in)  :: u_lab               (0:3,1:npsi,-nmu:nmu,i_s-2:i_e+2)
    real,     intent(in)  :: sin_row_theta
    real,     intent(in)  :: cos_row_theta
    integer,  intent(in)  :: direction

    real :: f_l(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real :: f_r(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real :: u_l_c    (0:3,1:npsi,-nmu:nmu)
    real :: u_r_c    (0:3,1:npsi,-nmu:nmu)

    real :: m_inv_c     (0:3,0:3)
    real :: lambda_inv_c(0:3,0:3)
    real :: mfull_inv_c (0:3,0:3)

    real :: left, center, right
    real :: slope_left, slope_right, slope_i

    real :: ir_if2(i_s-2:i_e+1)
    real :: r2    (i_s-2:i_e+2)
    real :: idx   (i_s-2:i_e+2)

    real :: dx    (i_s-1:i_e+2)

    integer :: a
    integer :: i, j, k
    integer :: l

    logical :: omp_in_parallel

    ! Calculate 1/r^2
    ! The central zone should not be divided by zero
    ! (doesn't matter because the area is zero anyway)
    if (direction == 1 .and. reconstruct_fr2) then
      ir_if2(i_s:i_e) = 1.0 / x_if(i_s:i_e)**2
      if (i_s == 1) then
        ir_if2(0)   = ir_if2(1) ! value at center wil lbe incorrect, but area is zero anyway
        ir_if2(-1)  = ir_if2(1)
      endif
      if (i_e == n) then
        ir_if2(n+1) = 1.0 / ( x_if(n) + (x_if(n) - x_if(n-1)) * (x_if(n) - x_if(n-1)) / (x_if(n-1) - x_if(n-2)) )**2
      endif
      r2 = x**2
    endif

    do a = i_s-1, i_e+2
      dx(a) = x(a) - x(a-1)
      idx(a) = 1.0 / dx(a)
    enddo

    !$omp parallel if (n > 1 .and. .not. omp_in_parallel()) &
    !$omp default(none) &
    !$omp private(a, i, j, k, l) &
    !$omp private(u_l_c, u_r_c) &
    !$omp private(f_l, f_r) &
    !$omp private(slope_left, slope_right, slope_i) &
    !$omp private(center, left, right) &
    !$omp private(m_inv_c) &
    !$omp private(lambda_inv_c, mfull_inv_c) &
    !$omp shared(n, i_s, i_e, direction, sin_row_theta, cos_row_theta) &
    !$omp shared(alpha, beta, phiconf) &
    !$omp shared(vfluid, x, x_if) &
    !$omp shared(u_lab) &
    !$omp shared(f, f_recon_l, f_recon_r) &
    !$omp shared(slope_l, slope_r) &
    !$omp shared(ir_if2, r2, idx, dx, reconstruct_fr2, slope_limiter)
    !$omp do
    do a = i_s-1, i_e+1
      ! Boost to comoving frame
      m_inv_c = m_inv_transform(alpha(a), beta(:,a), phiconf(a))
      lambda_inv_c = lambda_transform(-vfluid(:,a))
      mfull_inv_c = matmul(lambda_inv_c, m_inv_c)
      ! Bypass kick and approximate by going straight to next cell
      call boost(mfull_inv_c, u_lab(:,:,:,a-1), u_l_c)
      call boost(mfull_inv_c, u_lab(:,:,:,a+1), u_r_c)

      ! Transform left and right to center coordinate
      f_l = transform_frame(f(:,:,:,:,a-1), u_l_c,  dx(a),   sin_row_theta, cos_row_theta, direction)
      f_r = transform_frame(f(:,:,:,:,a+1), u_r_c, -dx(a+1), sin_row_theta, cos_row_theta, direction)

      do i = 1, neps
        do j = -nmu, nmu
          do k = 1, npsi
            do l = 1, nflav

              if (direction == 1 .and. reconstruct_fr2) then ! reconstruct fr^2
                center =   f(l,i,k,j,a) * r2(a)
                left   = f_l(l,i,k,j)   * r2(a-1)
                right  = f_r(l,i,k,j)   * r2(a+1)
              else
                center =   f(l,i,k,j,a)
                left   = f_l(l,i,k,j)
                right  = f_r(l,i,k,j)
              endif

              slope_left = (center - left) * idx(a)
              slope_right = (right - center) * idx(a+1)

              if (slope_limiter == 1) then

                ! Minmod limiter
                if (slope_left*slope_right > 0.0) then
                  if (abs(slope_left) < abs(slope_right)) then
                    slope_i = slope_left
                  else
                    slope_i = slope_right
                  endif
                else
                  slope_i = 0.0
                endif

              elseif (slope_limiter == 2) then

                ! MC limiter
                if (slope_left*slope_right > 0.0) then
                  slope_i = sign(1.0, slope_left) * &
                  min( abs(0.5*(slope_left+slope_right)), min(2.0*abs(slope_left), 2.0*abs(slope_right)) )
                else
                  slope_i = 0.0
                endif

              else ! no reconstruction

                slope_i = 0.0

              endif

              f_recon_l(l,i,k,j,a) = center + (x_if(a-1) - x(a))*slope_i
              f_recon_r(l,i,k,j,a) = center + (x_if(a  ) - x(a))*slope_i

              slope_l(l,i,k,j,a) = slope_i
              slope_r(l,i,k,j,a) = slope_i

              if (direction == 1 .and. reconstruct_fr2) then ! divide through by r^2
                slope_l(l,i,k,j,a) = slope_l(l,i,k,j,a) * ir_if2(a-1)
                slope_r(l,i,k,j,a) = slope_r(l,i,k,j,a) * ir_if2(a)
                if (a > 0  ) f_recon_l(l,i,k,j,a) = f_recon_l(l,i,k,j,a) * ir_if2(a-1)
                if (a < n+1) f_recon_r(l,i,k,j,a) = f_recon_r(l,i,k,j,a) * ir_if2(a)
              endif

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp enddo
    !$omp end parallel

  end subroutine reconstruct_linear


end module rice_reconstruct
