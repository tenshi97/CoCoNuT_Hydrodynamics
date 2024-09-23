module rice_utils
  implicit none
contains

  function determinant(matrix) result(det)
    real, intent(in) :: matrix(3,3)

    real :: det

    det = matrix(1,1) * (matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
        - matrix(1,2) * (matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
        + matrix(1,3) * (matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))

  end function determinant

  function cross_product(a, b) result(c)
    real, intent(in) :: a(3)
    real, intent(in) :: b(3)

    real :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product

  function magnitude(a) result(mag)
    real, intent(in) :: a(3)

    real :: mag

    mag = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
  end function magnitude

  function cartesian_coordinate(r, theta, phi) result(v)
    real, intent(in) :: r
    real, intent(in) :: theta
    real, intent(in) :: phi

    real :: v(3)

    v(1) = r * sin(theta) * cos(phi)
    v(2) = r * sin(theta) * sin(phi)
    v(3) = r * cos(theta)

  end function cartesian_coordinate

end module rice_utils
