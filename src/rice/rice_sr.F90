module rice_sr
  implicit none

contains

  function lambda_transform(v) result(lambda)
    use rice_config, only: clight
    real, intent(in)  :: v(3)
    real              :: lambda(0:3,0:3)

    real :: gamma
    real :: beta(3)
    real :: beta2
    real :: betai2

    beta(1) = v(1) / clight
    beta(2) = v(2) / clight
    beta(3) = v(3) / clight

    beta2 = dot_product(beta, beta)

    ! to stop divide-by-zero when velocity is zero
    if (abs(beta2) < epsilon(beta2)) then
      betai2 = 1.0
      gamma = 1.0
    else
      betai2 = 1.0 / beta2
      gamma = 1.0 / sqrt(1.0 - beta2)
    endif

    lambda(0,0) = gamma
    lambda(1,0) = gamma * beta(1)
    lambda(2,0) = gamma * beta(2)
    lambda(3,0) = gamma * beta(3)

    lambda(0,1) = lambda(1,0)
    lambda(1,1) = 1.0 + (gamma-1) * beta(1) * beta(1) * betai2
    lambda(2,1) =       (gamma-1) * beta(2) * beta(1) * betai2
    lambda(3,1) =       (gamma-1) * beta(3) * beta(1) * betai2

    lambda(0,2) = lambda(2,0)
    lambda(1,2) = lambda(2,1)
    lambda(2,2) = 1.0 + (gamma-1) * beta(2) * beta(2) * betai2
    lambda(3,2) =       (gamma-1) * beta(3) * beta(2) * betai2

    lambda(0,3) = lambda(3,0)
    lambda(1,3) = lambda(3,1)
    lambda(2,3) = lambda(3,2)
    lambda(3,3) = 1.0 + (gamma-1) * beta(3) * beta(3) * betai2

  end function lambda_transform

end module rice_sr
