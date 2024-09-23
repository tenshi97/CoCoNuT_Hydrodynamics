module rice_upwind
 implicit none

contains

  subroutine upwind_state(f_upw, f_l, f_r, u)
    use rice_config,       only: neps, nmu, npsi, nflav

    real,     intent(out) :: f_upw(1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: f_l  (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: f_r  (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real,     intent(in)  :: u    (1:nflav,1:neps,1:npsi,-nmu:nmu)

    integer :: i, j, k, l

    ! Select the left or right state based on upwinding using the velocity
    ! The velocity will be in the same direction as 'direction' variable
    do j = -nmu, nmu
      do k = 1, npsi
        do i = 1, neps
          do l = 1, nflav
            if (u(l,i,k,j) >= 0.0) then ! If velocity is zero, take inner cell
              f_upw(l,i,k,j) = f_l(l,i,k,j)
            else
              f_upw(l,i,k,j) = f_r(l,i,k,j)
            endif
          enddo
        enddo
      enddo
    enddo

  end subroutine upwind_state

end module rice_upwind
