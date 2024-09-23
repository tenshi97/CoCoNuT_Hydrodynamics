module rice_gausslobatto
  implicit none

contains

  subroutine legendre(n, x, pn, p1n, p2n)
  integer,  intent(in)    :: n
  real,     intent(in)    :: x
  real,     intent(out)   :: pn   ! legendre polynomial
  real,     intent(out)   :: p1n  ! first derivative
  real,     intent(out)   :: p2n  ! second derivative

  integer :: pmax

  real, allocatable :: p(:)
  real, allocatable :: p1(:)
  real, allocatable :: p2(:)

  integer :: j
  real    :: jr

  pmax = max(n ,3)

  allocate(p (0:pmax))
  allocate(p1(0:pmax))
  allocate(p2(0:pmax))

  p(0) = 0.0
  p(1) = 1.0
  p(2) = 0.5 * (3.0*x**2 - 1.0)
  p(3) = 0.5 * (5.0*x**3 - 3.0*x)

  p1(0) = 0.0
  p1(1) = 1.0
  p1(2) = 3.0 * x
  p1(3) = 0.5 * (15*x**2 - 3.0)

  p2(0) = 0.0
  p2(1) = 0.0
  p2(2) = 3.0
  p2(3) = 15.0 * x

  do j = 4, n
    jr = real(j)
    p(j)   = ((2.0*jr-1.0)*x*p(j-1) - (jr-1.0)*p(j-2)) / jr
    p1(j)  = jr*(x*p(j) - p(j-1)) / (x**2-1)
    p2(j)  = 4.0*x*p2(j-1) - (4.0*x**2+2.0)*p2(j-2) + 4.0*x*p2(j-3) - p2(j-4) + 2.0*p(j-2) + p1(j-1) - 2.0*x*p1(j-2) + p1(j-3)
  enddo

  pn  = p(n)
  p1n = p1(n)
  p2n = p2(n)

end subroutine legendre

subroutine collocation_points(nx, x, wt)
  integer,  intent(in)  :: nx
  real,     intent(out) :: x(1:nx)
  real,     intent(out) :: wt(1:nx)

  real :: p
  real :: p1
  real :: p2

  real :: err

  integer :: l
  integer :: k

  ! Edge x and weights
  x(1)   = -1.0
  x(nx) = 1.0
  wt(1)   = 2.0 / real(nx*(nx-1))
  wt(nx) = 2.0 / real(nx*(nx-1))

  ! Initialise guesses - this is always the RHS value as l increases
  x(2:nx-1) = 1.0

  ! Repeat for each lower order polynomial
  do l = 2, nx-1
    ! Initial guess as the half way point between the two roots of previous order
    x(2:l) = 0.5 * (x(1:l-1) + x(2:l))

    ! For each root
    do k = 2, l
      ! Newton-Raphson to find x(k)
      err = 1.0
      do while(abs(err) > epsilon(err))
        call legendre(l, x(k), p, p1, p2)
        err = p1/p2
        x(k) = x(k) - 0.5*err
      enddo
    enddo
  enddo

  ! Run once more to calculate weight
  do k = 2, nx-1
    call legendre(nx-1, x(k), p, p1, p2)
    wt(k) = 2.0 / (real(nx*(nx-1)) * p**2)
  enddo

end subroutine collocation_points

end module rice_gausslobatto
