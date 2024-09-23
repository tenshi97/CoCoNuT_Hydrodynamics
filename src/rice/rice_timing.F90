module rice_timing
  implicit none

  real :: tic
  real :: toc

contains

  subroutine start_timing
    use omp_lib
    tic = omp_get_wtime()
  end subroutine start_timing

  subroutine stop_timing
    use omp_lib
    toc = omp_get_wtime()
  end subroutine stop_timing

  subroutine print_timing(range, totsteps)
    integer, intent(in) :: range
    integer, intent(in) :: totsteps

    real :: time, time_per_step, time_total

    time = toc - tic
    time_per_step = time/real(range)
    time_total = time_per_step * real(totsteps)

    print*, "(timing) time per step:", time_per_step, "s"
    ! Double predicted time to account for slowdown
    print*, "(timing) predicted run time:", 2*int(time_total/60.0), "min"

  end subroutine print_timing

end module rice_timing
