module rice_ic
  implicit none

contains

  subroutine setup_ic()
    ! use rice_config, only: npsi, nmu, neps
    use rice_config, only: neps, nmu, npsi, nflav, rmin
    use rice_grid,   only: r, theta, f, f_eq, kappa_a, kappa_s, vfluid, time, jump, &
                           alpha, beta, phiconf, &
                           a_s, a_e, b_s, b_e, c_s, c_e
    use rice_output, only: read_output

    real :: height, radius

    integer :: a, b, c, i, j, k, l

    time = 0.0

    f = 0.0
    f_eq = 0.0
    kappa_a = 0.0
    kappa_s = 0.0
    vfluid = 0.0

    alpha   = 1.0
    beta    = 0.0
    phiconf = 1.0

    jump = .false.

    height = 0.0
    radius = 0.0

    ! Dummy check to supress warnings
    if (r(2) > r(1)) continue
    if (rmin > 0.0) continue
    if (theta(1) > 0.0) continue


    ! Empty loop to supress warnings
    do a = a_s, a_e
      do b = b_s, b_e
        do c = c_s, c_e
          do j = -nmu, nmu
            do k = 1, npsi
              do i = 1, neps
                do l = 1, nflav
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! TEST 1
    ! print*, '(ic) Zero initial conditions'
    ! f = 0.0
    ! END TEST

    ! TEST 2
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e
          ! do j = -nmu, nmu
            ! do k = 1, npsi
              ! do i = 1, neps
                ! do l = 1, nflav
                  ! if (j == -nmu) then
                    ! if (0.8 < r(a) .and. r(a) < 1.2) then
                      ! f(l,i,k,j,c,b,a) = exp(-(r(a) - 1.0)**2 / 0.05**2)
                    ! endif
                  ! endif
                ! enddo
              ! enddo
            ! enddo
          ! enddo
        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 3
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e
          ! if (1.2 < r(a) .and. r(a) < 2.0) then
            ! vfluid(1,c,b,a) = -0.3
          ! endif
          ! do j = -nmu, nmu
            ! do k = 1, npsi
              ! do i = 1, neps
                ! do l = 1, nflav
                  ! if (j == nmu .and. i == neps/2) then
                    ! if (0.5 < r(a) .and. r(a) < 1.0) then
                      ! f(l,i,k,j,c,b,a) = exp(-(r(a) - 1.0)**2 / 0.05**2)
                    ! endif
                  ! endif
                ! enddo
              ! enddo
            ! enddo
          ! enddo
        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 4
    ! print*, '(ic) Radiating sphere initial conditions'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! if (r(a) < 1.0) then
            ! f_eq(1,1,c,b,a) = 1.0
            ! kappa_a(1,1,c,b,a) = 100.0
          ! endif

        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_04.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 5
    ! print*, '(ic) Redshift initial conditions'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e
          ! do i = 1, neps

            ! if (r(a) < 1.0 .and. i == neps / 2) then
              ! f_eq(1,i,c,b,a) = 1.0
              ! kappa_a(1,i,c,b,a) = 100.0
            ! endif

          ! enddo

          ! if (2.0 < r(a) .and. r(a) < 4.0) then
            ! vfluid(1,c,b,a) = - 0.3 * 2.0 / r(a)
            ! vfluid(1,c,b,a) = - 0.3 * exp(-(r(a) - 3.0)**2 / 0.4**2)
          ! endif

        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 6
    ! print*, '(ic) Pre-SN initial conditions'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! vfluid(1,c,b,a) = - 0.3 * exp(-(log10(r(a)) - 8.0)**2 / 1.0**2)

          ! do i = 1, neps

            ! if (r(a) < 1.0e8) then
              ! f_eq(1,i,c,b,a) = 1.0
              ! kappa_a(1,i,c,b,a) = 1.e-10 * real(i)*real(i)
            ! endif

          ! enddo

        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 7
    ! print*, '(ic) 2D radiating sphere'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e
          ! if (r(a) <= 1.0) then
            ! f_eq   (1,1,c,b,a) = 1.0
            ! kappa_a(1,1,c,b,a) = 1000.0
          ! endif
        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_07.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 8
    ! print*, '(ic) GR radiating sphere initial conditions varying with alpha and phiconf'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! if (r(a) < 1.0) then
            ! f_eq(1,neps,c,b,a) = 1.0
            ! kappa_a(1,:,c,b,a) = 100.0 ! Apply to all energy groups to prevent core build-up
          ! endif

          ! alpha(c,b,a) = 0.5
          ! if (r(a) >= 1.2) then
            ! alpha(c,b,a) = alpha(c,b,a) &
                         ! + (1.0 - alpha(c,b,a)) * (1.0 - 1.2/r(a)) &
                         ! * (1.0 - exp(-(1.2-r(a))**2/0.1))
          ! endif

          ! phiconf(c,b,a) = 1.0 / sqrt(alpha(c,b,a))

        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_08.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 9
    ! print*, '(ic) GR radiating sphere initial conditions varying just alpha'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! if (r(a) < 1.0) then
            ! f_eq(1,neps,c,b,a) = 1.0
            ! kappa_a(1,:,c,b,a) = 100.0 ! Apply to all energy groups to prevent core build-up
          ! endif

          ! alpha(c,b,a) = 0.5
          ! if (r(a) >= 1.2) then
            ! alpha(c,b,a) = alpha(c,b,a) &
                         ! + (1.0 - alpha(c,b,a)) * (1.0 - 1.2/r(a)) &
                         ! * (1.0 - exp(-(1.2-r(a))**2/0.1))
          ! endif

          ! phiconf(c,b,a) = 1.0 / sqrt(alpha(c,b,a))
          ! phiconf(c,b,a) = 1.0

        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_09.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 10
    ! print*, '(ic) GR radiating sphere initial conditions varying just phiconf'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! if (r(a) < 1.0) then
            ! f_eq(1,neps,c,b,a) = 1.0
            ! kappa_a(1,:,c,b,a) = 100.0 ! Apply to all energy groups to prevent core build-up
          ! endif

          ! alpha(c,b,a) = 0.5
          ! if (r(a) >= 1.2) then
            ! alpha(c,b,a) = alpha(c,b,a) &
                         ! + (1.0 - alpha(c,b,a)) * (1.0 - 1.2/r(a)) &
                         ! * (1.0 - exp(-(1.2-r(a))**2/0.1))
          ! endif

        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_10.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 11
    ! print*, '(ic) 2D conservation test'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e
          ! if (r(a) <= 1.0) then
            ! f(:,:,:,:,c,b,a) = 1.0
          ! endif
        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 12
    ! print*, '(ic) 2D radiating disk'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! radius = r(a) * sin(theta(b))
          ! height = r(a) * cos(theta(b))
          ! if (height <= 0.0) then
            ! if (r(a) <= 1.0) then
              ! f_eq   (1,1,c,b,a) = 1.0
              ! kappa_a(1,1,c,b,a) = 1.e6
            ! else
              ! f_eq   (1,1,c,b,a) = 0.0
              ! kappa_a(1,1,c,b,a) = 1.e6
            ! endif
          ! endif

        ! enddo
      ! enddo
    ! enddo
    ! call read_output('tests/ic/test_12.hdf5') ! Read a stationary run to save time
    ! END TEST

    ! TEST 13
    ! print*, '(ic) Performance test (same as test 8, without reading ICs from file)'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! if (r(a) < 1.0) then
            ! f_eq(1,neps,c,b,a) = 1.0
            ! kappa_a(1,:,c,b,a) = 100.0 ! Apply to all energy groups to prevent core build-up
          ! endif

          ! alpha(c,b,a) = 0.5
          ! if (r(a) >= 1.2) then
            ! alpha(c,b,a) = alpha(c,b,a) &
                         ! + (1.0 - alpha(c,b,a)) * (1.0 - 1.2/r(a)) &
                         ! * (1.0 - exp(-(1.2-r(a))**2/0.1))
          ! endif

          ! phiconf(c,b,a) = 1.0 / sqrt(alpha(c,b,a))

        ! enddo
      ! enddo
    ! enddo
    ! END TEST

    ! TEST 14
    ! print*, '(ic) 2D MPI identity test'
    ! do a = a_s, a_e
      ! do b = b_s, b_e
        ! do c = c_s, c_e

          ! do j = -nmu, nmu
            ! do k = 1, npsi
              ! do i = 1, neps
                ! do l = 1, nflav
                  ! f(l,i,k,j,c,b,a) = a + b + c + i + j + k + l
                ! enddo
              ! enddo
            ! enddo
          ! enddo

        ! enddo
      ! enddo
    ! enddo
    ! END TEST

  end subroutine setup_ic

end module rice_ic
