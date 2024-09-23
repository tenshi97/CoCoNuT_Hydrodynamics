module rice_boost
  implicit none

contains

  subroutine boost(m, u, uboost)
    use rice_config, only: nmu, npsi
    real, intent(in)  :: m (0:3,0:3)
    real, intent(in)  :: u     (0:3,1:npsi,-nmu:nmu)
    real, intent(out) :: uboost(0:3,1:npsi,-nmu:nmu)

    integer :: j, k

    do j = -nmu, nmu
      do k = 1, npsi
        uboost(:,k,j) = matmul(m, u(:,k,j))
      enddo
    enddo

  end subroutine boost

  subroutine boost_split(m, u, u_adv, u_trans)
    use rice_config, only: nmu, npsi
    real, intent(in)  :: m  (0:3,0:3)
    real, intent(in)  :: u      (0:3,1:npsi,-nmu:nmu)
    real, intent(out) :: u_adv  (0:3,1:npsi,-nmu:nmu)
    real, intent(out) :: u_trans(0:3,1:npsi,-nmu:nmu)

    integer :: j, k

    do j = -nmu, nmu
      do k = 1, npsi
        ! Fluid advection component
        u_adv(0,k,j) = m(0,0)*u(0,k,j)
        u_adv(1,k,j) = m(1,0)*u(0,k,j)
        u_adv(2,k,j) = m(2,0)*u(0,k,j)
        u_adv(3,k,j) = m(3,0)*u(0,k,j)

        ! Transport component
        u_trans(0,k,j) = m(0,1)*u(1,k,j) + m(0,2)*u(2,k,j) + m(0,3)*u(3,k,j)
        u_trans(1,k,j) = m(1,1)*u(1,k,j) + m(1,2)*u(2,k,j) + m(1,3)*u(3,k,j)
        u_trans(2,k,j) = m(2,1)*u(1,k,j) + m(2,2)*u(2,k,j) + m(2,3)*u(3,k,j)
        u_trans(3,k,j) = m(3,1)*u(1,k,j) + m(3,2)*u(2,k,j) + m(3,3)*u(3,k,j)

      enddo
    enddo

  end subroutine boost_split

end module rice_boost
