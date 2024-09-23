module rice_gr
  implicit none

contains

  function metric(alpha, beta, phiconf) result(g)
    real, intent(in)  :: alpha
    real, intent(in)  :: beta(3)
    real, intent(in)  :: phiconf

    real :: g(0:3,0:3)

    real :: beta_down(3)
    real :: phiconf4

    phiconf4 = phiconf**4
    beta_down(:) = beta*phiconf4

    g(0,0) = -alpha**2 + beta_down(1)*beta(1) + beta_down(2)*beta(2) + beta_down(3)*beta(3)
    g(0,1) = beta_down(1)
    g(0,2) = beta_down(2)
    g(0,3) = beta_down(3)

    g(1,0) = beta_down(1)
    g(1,1) = phiconf4
    g(1,2) = 0.0
    g(1,3) = 0.0

    g(2,0) = beta_down(2)
    g(2,1) = 0.0
    g(2,2) = phiconf4
    g(2,3) = 0.0

    g(3,0) = beta_down(3)
    g(3,1) = 0.0
    g(3,2) = 0.0
    g(3,3) = phiconf4

  end function metric

  function metric_inv(alpha, beta, phiconf) result(g_inv)
    real, intent(in)  :: alpha
    real, intent(in)  :: beta(3)
    real, intent(in)  :: phiconf

    real :: g_inv(0:3,0:3)

    real :: alphai2
    real :: phiconfi4

    alphai2 = 1.0 / alpha**2
    phiconfi4 = 1.0 / phiconf**4

    g_inv(0,0) = - alphai2
    g_inv(0,1) =   alphai2 * beta(1)
    g_inv(0,2) =   alphai2 * beta(2)
    g_inv(0,3) =   alphai2 * beta(3)

    g_inv(1,0) =   alphai2 * beta(1)
    g_inv(1,1) =   phiconfi4 - alphai2*beta(1)*beta(1)
    g_inv(1,2) =             - alphai2*beta(1)*beta(2)
    g_inv(1,3) =             - alphai2*beta(1)*beta(3)

    g_inv(2,0) =   alphai2 * beta(2)
    g_inv(2,1) =             - alphai2*beta(2)*beta(1)
    g_inv(2,2) =   phiconfi4 - alphai2*beta(2)*beta(2)
    g_inv(2,3) =             - alphai2*beta(2)*beta(3)

    g_inv(3,0) =   alphai2 * beta(3)
    g_inv(3,1) =             - alphai2*beta(3)*beta(1)
    g_inv(3,2) =             - alphai2*beta(3)*beta(2)
    g_inv(3,3) =   phiconfi4 - alphai2*beta(3)*beta(3)

  end function metric_inv

  function m_transform(alpha, beta, phiconf) result(m)
    real, intent(in)  :: alpha
    real, intent(in)  :: beta(3)
    real, intent(in)  :: phiconf

    real :: m(0:3,0:3)

    real :: alphai
    real :: phiconfi2

    alphai = 1.0 / alpha
    phiconfi2 = 1.0 / phiconf**2

    m(0,0) = alphai
    m(0,1) = 0.0
    m(0,2) = 0.0
    m(0,3) = 0.0

    m(1,0) = alphai * beta(1)
    m(1,1) = phiconfi2
    m(1,2) = 0.0
    m(1,3) = 0.0

    m(2,0) = alphai * beta(2)
    m(2,1) = 0.0
    m(2,2) = phiconfi2
    m(2,3) = 0.0

    m(3,0) = alphai * beta(3)
    m(3,1) = 0.0
    m(3,2) = 0.0
    m(3,3) = phiconfi2

  end function m_transform

  function m_inv_transform(alpha, beta, phiconf) result(m_inv)
    real, intent(in)  :: alpha
    real, intent(in)  :: beta(3)
    real, intent(in)  :: phiconf

    real :: m_inv(0:3,0:3)

    real :: phiconf2

    phiconf2 = phiconf**2

    m_inv(0,0) = alpha
    m_inv(0,1) = phiconf2 * beta(1)
    m_inv(0,2) = phiconf2 * beta(2)
    m_inv(0,3) = phiconf2 * beta(3)

    m_inv(1,0) = 0.0
    m_inv(1,1) = phiconf2
    m_inv(1,2) = 0.0
    m_inv(1,3) = 0.0

    m_inv(2,0) = 0.0
    m_inv(2,1) = 0.0
    m_inv(2,2) = phiconf2
    m_inv(2,3) = 0.0

    m_inv(3,0) = 0.0
    m_inv(3,1) = 0.0
    m_inv(3,2) = 0.0
    m_inv(3,3) = phiconf2

  end function m_inv_transform

  function kick(g_inv_lr, dg, u_if, i) result(delta)
    real,    intent(in) :: g_inv_lr(0:3,0:3)
    real,    intent(in) :: dg      (0:3,0:3)
    real,    intent(in) :: u_if(0:3)
    integer, intent(in) :: i

    real :: delta
    real :: a, b, disc

    call kick_quadratic_terms(g_inv_lr, dg, u_if, i, a=a, b=b, disc=disc)

    ! If disc < 0, there are no real solutions, return maximal kick in i direction
    if (disc < 0.0 .and. i==1) then
      delta = -u_if(i)
      return
    endif

    ! Quadratic equation, factor of 2 taken out
    ! If velocity is positive, take positive solution
    if (u_if(i) >= 0.0 .or. i == 0) then
      delta = (-b + sqrt(disc)) / a
    ! If velocity is negative, take negative solution
    else
      delta = (-b - sqrt(disc)) / a
    endif

  end function kick

  function modified_kick(g_inv_lr, g_inv_if, u_if) result(k)
    ! only works for i=1
    real,    intent(in) :: g_inv_lr(0:3,0:3)
    real,    intent(in) :: g_inv_if(0:3,0:3)
    real,    intent(in) :: u_if(0:3)

    real :: k
    real :: a, b, c
    real :: guu_if

    guu_if = &
    + g_inv_if(0,0)*u_if(0)*u_if(0) &
    + 2.0*g_inv_if(1,0)*u_if(1)*u_if(0) + g_inv_if(1,1)*u_if(1)*u_if(1) &
    + 2.0*g_inv_if(2,0)*u_if(2)*u_if(0) + 2.0*g_inv_if(2,1)*u_if(2)*u_if(1) + g_inv_if(2,2)*u_if(2)*u_if(2) &
    + 2.0*g_inv_if(3,0)*u_if(3)*u_if(0) + 2.0*g_inv_if(3,1)*u_if(3)*u_if(1) + 2.0*g_inv_if(3,2)*u_if(3)*u_if(2) &
                                                                            + g_inv_if(3,3)*u_if(3)*u_if(3)

    a =     g_inv_lr(2,2)*u_if(2)*u_if(2) + 2.0*g_inv_lr(2,3)*u_if(2)*u_if(3) + g_inv_lr(3,3)*u_if(3)*u_if(3)
    b = 2.0*g_inv_lr(0,2)*u_if(0)*u_if(2) + 2.0*g_inv_lr(0,3)*u_if(0)*u_if(3)
    c =     g_inv_lr(0,0)*u_if(0)*u_if(0) - guu_if

    k = (-b + sqrt(b**2 - 4.0*a*c)) / (2.0*a)

  end function modified_kick

  subroutine kick_quadratic_terms(g_inv_lr, dg, u_if, i, a, b, c, disc)
    real,           intent(in)  :: g_inv_lr(0:3,0:3)
    real,           intent(in)  :: dg      (0:3,0:3)
    real,           intent(in)  :: u_if(0:3)
    integer,        intent(in)  :: i
    real, optional, intent(out) :: a
    real, optional, intent(out) :: b
    real, optional, intent(out) :: c
    real, optional, intent(out) :: disc

    real :: a_, b_, c_

    if (present(a) .or. present(disc)) then
      a_ = g_inv_lr(i,i)
    endif

    if (present(b) .or. present(disc)) then
      b_ = u_if(0)*g_inv_lr(0,i) + u_if(1)*g_inv_lr(1,i) + u_if(2)*g_inv_lr(2,i) + u_if(3)*g_inv_lr(3,i)
    endif

    if (present(c) .or. present(disc)) then
      c_ = dg(0,0)*u_if(0)*u_if(0) &
         + 2.0*dg(1,0)*u_if(1)*u_if(0) + dg(1,1)*u_if(1)*u_if(1) &
         + 2.0*dg(2,0)*u_if(2)*u_if(0) + 2.0*dg(2,1)*u_if(2)*u_if(1) + dg(2,2)*u_if(2)*u_if(2) &
         + 2.0*dg(3,0)*u_if(3)*u_if(0) + 2.0*dg(3,1)*u_if(3)*u_if(1) + 2.0*dg(3,2)*u_if(3)*u_if(2) + dg(3,3)*u_if(3)*u_if(3)
    endif

    if (present(a)) then
      a = a_
    endif

    if (present(b)) then
      b = b_
    endif

    if (present(c)) then
      c = c_
    endif

    if (present(disc)) then
      disc = b_**2 - a_*c_
    endif

  end subroutine kick_quadratic_terms

  subroutine ricci_rotation(g_src, g_inv_src, g_inv_dst, u_src, u_dst, direction)
    use rice_config, only: nmu, npsi
    use rice_boost,  only: boost
    real,    intent(in)  :: g_src    (0:3,0:3)
    real,    intent(in)  :: g_inv_src(0:3,0:3)
    real,    intent(in)  :: g_inv_dst(0:3,0:3)
    real,    intent(in)  :: u_src(0:3,1:npsi,-nmu:nmu)
    real,    intent(out) :: u_dst(0:3,1:npsi,-nmu:nmu)
    integer, intent(in)  :: direction

    real :: u_down(0:3,1:npsi,-nmu:nmu)
    real :: dg(0:3,0:3)

    integer :: j, k

    real :: kfactor, u_new

    ! No rotation in lateral directions
    if (direction == 2 .or. direction == 3) then
      u_dst = u_src
      return
    endif

    call boost(g_src, u_src, u_down)
    dg = g_inv_dst - g_inv_src

    do j = -nmu, nmu
      do k = 1, npsi

        ! Get new value after kick
        u_new = u_down(direction,k,j) + kick(g_inv_src, dg, u_down(:,k,j), direction)

        ! If u_new = 0, then no conservative solution exists
        if (direction == 1 .and. .not. (u_new >= 0.0 .eqv. u_down(direction,k,j) >= 0.0)) then
          ! Do energy (but not momentum) conservative kick
          kfactor = modified_kick(g_inv_dst, g_inv_src, u_down(:,k,j))
          u_down(1,k,j) = 0.0
          u_down(2,k,j) = kfactor * u_down(2,k,j)
          u_down(3,k,j) = kfactor * u_down(3,k,j)
        else
          ! Do fully conservative kick
          ! For all directions > 1, disc should >= 0
          u_down(direction,k,j) = u_new
        endif

      enddo
    enddo

    call boost(g_inv_dst, u_down, u_dst)

  end subroutine ricci_rotation

  subroutine gr_flux_correction(flux, phiconf, phiconf_if, alpha, alpha_if)
    use rice_config, only: nflav, neps, npsi, nmu
    real, intent(inout) :: flux(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(in)    :: phiconf
    real, intent(in)    :: phiconf_if
    real, intent(in)    :: alpha
    real, intent(in)    :: alpha_if

    flux = flux * (alpha_if/alpha) * (phiconf_if/phiconf)**6

  end subroutine gr_flux_correction

  subroutine gr_flux_u0_divide(flux,  u)
    use rice_config, only: nflav, neps, npsi, nmu
    real, intent(inout) :: flux(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(in)    :: u(0:3,1:npsi,-nmu:nmu)

    integer :: j, k

    do j = -nmu, nmu
      do k = 1, npsi
        flux(:,:,k,j) = flux(:,:,k,j) / u(0,k,j)
      enddo
    enddo

  end subroutine gr_flux_u0_divide

  subroutine critical_u_down(g_inv_dst, dg, u_down, u_crit, tir)
    real,    intent(in)  :: g_inv_dst(0:3,0:3)
    real,    intent(in)  :: dg       (0:3,0:3)
    real,    intent(in)  :: u_down   (0:3)
    real,    intent(out) :: u_crit   (0:3)
    logical, intent(out) :: tir

    integer :: i

    real :: p, p2
    real :: r, r2, rnew, ri
    real :: x
    real :: disc, ddisc, d2disc
    real :: dx

    r2 = u_down(2)**2 + u_down(3)**2

    ! If mu = -1 or 1
    if (r2 < epsilon(r2)) then
      tir = .false.
      return
    endif

    ! Check if u can cross conservatively
    call kick_quadratic_terms(g_inv_dst, dg, u_down, 1, disc=disc)
    if (disc >= 0.0)  then
      tir = .false.
      return
    endif

    ! If disc < 0:
    tir = .true.
    u_crit = u_down

    r = sqrt(r2)
    ri = 1.0 / r
    p2 = u_down(1)**2 + r2
    p = sqrt(p2)

    ! Starting guess
    x = 0.01

    do i = 1, 100
      ! Get new disc
      call ddiscdx(x, p, g_inv_dst, dg, u_crit(0), disc, ddisc, d2disc)

      ! If within tolerance, exit
      if (disc > 0.0) then
        if (disc < 2.0*epsilon(disc) .or. abs(ddisc) < epsilon(ddisc)) exit
      endif

      ! Newtopn-Raphson
      ! dx = (disc - epsilon(disc)) / ddisc

      ! Use 2nd derivative
      dx = (-ddisc + sqrt(ddisc**2 + 2.0*d2disc*(disc-epsilon(disc)))) / d2disc

      x = x - dx

    enddo

    ! Get new vector
    rnew = p * sqrt(1.0 - x**2) * ri
    u_crit(1) = sign(x * p, u_down(1))
    u_crit(2) = u_down(2) * rnew
    u_crit(3) = u_down(3) * rnew

  end subroutine critical_u_down

  subroutine ddiscdx(x, p, g, dg, u0, disc, ddisc, d2disc)
    real, intent(in) :: x
    real, intent(in) :: p
    real, intent(in) :: g (0:3,0:3)
    real, intent(in) :: dg(0:3,0:3)
    real, intent(in) :: u0
    real, intent(out) :: disc
    real, intent(out) :: ddisc
    real, intent(out) :: d2disc

    real :: x2
    real :: mx, sx, sxi, sx3
    real :: px
    real :: term

    x2 = x**2
    mx = 1.0 - x2
    sx = sqrt(mx)
    sx3 = sx*mx
    sxi = 1.0 / sx
    px = p*x

    disc = (g(0,1)*u0 + g(1,1)*px + g(1,2)*p*sx)**2 - g(1,1) * (dg(0,0)*u0**2 &
        + dg(2,2)*p**2*sx**2 + p*(2.0*dg(0,1)*u0*x + dg(1,1)*px*x + 2.0*(dg(0,2)*u0 + dg(1,2)*px)*sx))

    ddisc = 2.0 * p * ( &
                       dg(2,2)*g(1,1)*px &
                     - g(1,1)*(dg(0,1)*u0 + dg(1,1)*px - x*(dg(0,2)*u0 + dg(1,2)*px)*sxi + dg(1,2)*p*sxi) &
                     + (g(1,1) - g(1,2)*x*sxi) * (g(0,1)*u0 + g(1,1)*px + g(1,2)*p*sx) &
                     )

    term = x*(2.0*x2 - 3.0)

    d2disc = 2.0 * p * ( &
                      - dg(0,2)*g(1,1)*u0 + g(0,1)*g(1,2)*u0 &
                      + p * ( &
                              (dg(1,1) - dg(2,2))*g(1,1)*sx3 &
                            - g(1,1)**2*sx3 &
                            + g(1,2)**2*sx3 &
                            + dg(1,2)*g(1,1)*term &
                            -  g(1,1)*g(1,2)*term &
                            ) &
                      )

  end subroutine ddiscdx
end module rice_gr
