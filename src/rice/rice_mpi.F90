module rice_mpi
  use rice_grid, only: sbuf, rbuf, &
                       a_s, a_e, b_s, b_e, c_s, c_e, &
                       nbuf, ibuf, ibuf_send
#ifdef MPI_BT
#ifdef USE_MPIF_H
  implicit none
  include "mpif.h"
#else
  use mpi
  implicit none
#endif
#else
  implicit none
#endif

  integer :: id
  integer :: cart_comm
  integer :: mpi_dims(3)
  integer :: ncyc(3)

  interface fill_buffer
    module procedure fill_buffer_scalar, &
                      fill_buffer_1d, &
                      fill_buffer_2d, &
                      fill_buffer_3d, &
                      fill_buffer_4d
  end interface

  interface read_buffer
    module procedure read_buffer_scalar, &
                      read_buffer_1d, &
                      read_buffer_2d, &
                      read_buffer_3d, &
                      read_buffer_4d
  end interface

contains

  subroutine setup_mpi(init)
    use rice_config, only: mpi_nr, mpi_ntheta, mpi_nphi
#ifdef MPI_BT
    use rice_grid,   only: cart_coords, grid_2d
#endif

    logical, optional :: init

#ifndef MPI_BT
    if (mpi_nr > 1)     stop 'MPI not enabled, but r domain decomposition requested'
    if (mpi_ntheta > 1) stop 'MPI not enabled, but theta domain decomposition requested'
    if (mpi_nphi > 1)   stop 'MPI not enabled, but phi domain decomposition requested'

    id = 0

    ! Supress warning
    if (present(init)) continue

#else

    integer :: cart_dims(0:2)
    logical :: cart_periods(0:2)
    integer :: ierr

    if (present(init)) then
      if (init) call MPI_Init(ierr)
    endif

    call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)
    print*, '(mpi) rank', id

    mpi_dims(1) = mpi_nr
    mpi_dims(2) = mpi_ntheta
    mpi_dims(3) = mpi_nphi

    print*, '(mpi) domain decomposition:', mpi_dims

    cart_dims(0:2) = mpi_dims(1:3)

    ! Physical BCs are applied later
    cart_periods(0) = .false.
    cart_periods(1) = (.not. grid_2d)
    cart_periods(2) = .true.

    print*, '(mpi) creating cartesian communicator'
    call MPI_Cart_create(MPI_COMM_WORLD, 3, cart_dims, cart_periods, .false., cart_comm, ierr)
    call MPI_Cart_coords(cart_comm, id, 3, cart_coords, ierr)

#endif

  end subroutine setup_mpi

  subroutine finish_mpi
#ifdef MPI_BT
    integer :: ierr
    call MPI_Finalize(ierr)
#endif
  end subroutine finish_mpi

  subroutine calculate_ncyc
    use rice_config, only: mpi_nr, mpi_ntheta, mpi_nphi, nr, ntheta, nphi
    ! If there is only one cell-width on the task, then 2 MPI communication cycles
    ! are required to propgate ghost cells
    ncyc(:) = 1
    if (nr     == mpi_nr    ) ncyc(1) = 2
    if (ntheta == mpi_ntheta) ncyc(2) = 2
    if (nphi   == mpi_nphi  ) ncyc(3) = 2
    print*, '(mpi) number of cycles in each direction:', ncyc
  end subroutine calculate_ncyc

  subroutine exchange

    integer :: direction, side
    integer :: cyc
#ifdef MPI_BT
    integer :: tag
    integer :: src, dst, ierr
    integer :: mpistat(MPI_STATUS_SIZE)
#endif

    do direction = 1, 3
      do side = 0, 1
        do cyc = 1, ncyc(direction)
          call fill_send_buffer(side, direction)
#ifdef MPI_BT
          ! CoCoNuT only sets up the cartesian communicator for
          ! theta: cart direction 0
          ! phi  : cart direction 1
          !
          ! for r direction, call self_buffer_copy
          !
          ! disp = -1: send to the left, receive from the right
          ! disp = +1: send to the right, receive from the left
          !
          ! side = 0 fills from the left and unfills on the right
          ! so send this to the left (-1)
          ! side = 1 fills from the right and unfills on the left
          ! so send this to the right (+1)

          if (mpi_dims(direction) == 1) then
            ! If there is only 1 task in this direction, then do a direct memory copy
            call self_buffer_copy
          else
            call MPI_Cart_shift(cart_comm, direction-1, -1+side*2, src, dst, ierr)

            tag = direction * 10 + side
            call MPI_Sendrecv(sbuf(1:ibuf), ibuf, MPI_DOUBLE_PRECISION, dst, tag, &
                              rbuf(1:ibuf), ibuf, MPI_DOUBLE_PRECISION, src, tag, &
                              cart_comm, mpistat, ierr)
          endif

#else
          ! If not compiled with MPI, this routine simply exchanges
          ! the LHS and RHS physical ghosts because there is only 1 task

          ! call self_buffer_copy

          ! Instead of copying memory, use pointers. rbuf => sbuf
#endif
          call read_recv_buffer(mod(side+1, 2), direction)

        enddo
      enddo
    enddo

  end subroutine exchange

  subroutine self_buffer_copy

    ! Copy the send buffer directly into the receive buffer
    ! for serial runs
    rbuf(1:ibuf) = sbuf(1:ibuf)

  end subroutine self_buffer_copy

  subroutine fill_send_buffer(side, direction)
    use rice_config, only: neps, nmu, npsi, nflav
    use rice_grid, only: f_eq, kappa_a, kappa_s, vfluid, alpha, beta, phiconf, f, nbuf
    integer, intent(in) :: side
    integer, intent(in) :: direction

    ! Fills the send buffer with:
    ! 1,   2 for LHS
    ! n-1, n for RHS
    ibuf = 0
    call fill_buffer(f_eq,                     neps, nflav, side, direction)
    call fill_buffer(kappa_a,                  neps, nflav, side, direction)
    call fill_buffer(kappa_s,                  neps, nflav, side, direction)
    call fill_buffer(vfluid,                             3, side, direction)
    call fill_buffer(alpha,                                 side, direction)
    call fill_buffer(beta,                               3, side, direction)
    call fill_buffer(phiconf,                               side, direction)
    call fill_buffer(f,       (2*nmu+1), npsi, neps, nflav, side, direction)

    ! Check entire buffer has been filled
    if (ibuf > nbuf) stop 'MPI buffer write overflow'

    ibuf_send = ibuf

  end subroutine fill_send_buffer

  subroutine fill_buffer_scalar(x, side, direction)
    real,    intent(in)    :: x(c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .false., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      sbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1) = &
          RESHAPE(x(c_s_buf:c_e_buf,b_s_buf:b_e_buf,a),(/nn/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine fill_buffer_scalar

  subroutine fill_buffer_1d(x, n1, side, direction)
    integer, intent(in)    :: n1
    real,    intent(in)    :: x(n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .false., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      sbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1) = &
          RESHAPE(x(1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a),(/nn/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine fill_buffer_1d

  subroutine fill_buffer_2d(x, n1, n2, side, direction)
    integer, intent(in)    :: n1, n2
    real,    intent(in)    :: x(n2,n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .false., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      sbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1) = &
          RESHAPE(x(1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a),(/nn/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine fill_buffer_2d

  subroutine fill_buffer_3d(x, n1, n2, n3, side, direction)

    implicit none

    integer, intent(in)    :: n1, n2, n3
    real,    intent(in)    :: x(n3,n2,n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .false., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2 * n3

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      sbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1) = &
          RESHAPE(x(1:n3,1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a),(/nn/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn- 1


  end subroutine fill_buffer_3d

  subroutine fill_buffer_4d(x, n1, n2, n3, n4, side, direction)

    implicit none

    integer, intent(in)    :: n1, n2, n3, n4
    real,    intent(in)    :: x(n4,n3,n2,n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .false., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2 * n3 * n4

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      sbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1) = &
          RESHAPE(x(1:n4,1:n3,1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a),(/nn/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine fill_buffer_4d

  subroutine read_recv_buffer(side, direction)
    use rice_config, only: neps, nmu, npsi, nflav
    use rice_grid, only: f_eq, kappa_a, kappa_s, vfluid, alpha, beta, phiconf, f, nbuf
    integer, intent(in) :: side
    integer, intent(in) :: direction

    ! Reads the recv buffer into:
    ! -1,  0   for LHS
    ! n+1, n+2 for RHS

    ! Read must be performed in exactly the same order as fill
    ! otherwise the data will be read into the wrong arrays
    ibuf = 0
    call read_buffer(f_eq,                     neps, nflav, side, direction)
    call read_buffer(kappa_a,                  neps, nflav, side, direction)
    call read_buffer(kappa_s,                  neps, nflav, side, direction)
    call read_buffer(vfluid,                             3, side, direction)
    call read_buffer(alpha,                                 side, direction)
    call read_buffer(beta,                               3, side, direction)
    call read_buffer(phiconf,                               side, direction)
    call read_buffer(f,       (2*nmu+1), npsi, neps, nflav, side, direction)

    ! Sanity checks
    if (ibuf > nbuf) stop 'MPI buffer read overflow'
    if (ibuf /= ibuf_send) stop 'MPI read data is not the same length as send data'

  end subroutine read_recv_buffer

  subroutine read_buffer_scalar(x, side, direction)
    real,    intent(inout) :: x(c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .true., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      x(c_s_buf:c_e_buf,b_s_buf:b_e_buf,a) = &
          RESHAPE(rbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1), (/nc,nb/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine read_buffer_scalar

  subroutine read_buffer_1d(x, n1, side, direction)
    integer, intent(in)    :: n1
    real,    intent(inout) :: x(1:n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .true., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      x(1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a) = &
          RESHAPE(rbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1), (/n1,nc,nb/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine read_buffer_1d

  subroutine read_buffer_2d(x, n1, n2, side, direction)
    integer, intent(in)    :: n1
    integer, intent(in)    :: n2
    real,    intent(inout) :: x(1:n2,1:n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .true., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      x(1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a) = &
          RESHAPE(rbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1), (/n2,n1,nc,nb/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine read_buffer_2d

  subroutine read_buffer_3d(x, n1, n2, n3, side, direction)

    implicit none

    integer, intent(in)    :: n1
    integer, intent(in)    :: n2
    integer, intent(in)    :: n3
    real,    intent(inout) :: x(1:n3,1:n2,1:n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .true., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2 * n3

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      x(1:n3,1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a) = &
          RESHAPE(rbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1), (/n3,n2,n1,nc,nb/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine read_buffer_3d

  subroutine read_buffer_4d(x, n1, n2, n3, n4, side, direction)

    implicit none

    integer, intent(in)    :: n1
    integer, intent(in)    :: n2
    integer, intent(in)    :: n3
    integer, intent(in)    :: n4
    real,    intent(inout) :: x(1:n4,1:n3,1:n2,1:n1,c_s-2:c_e+2,b_s-2:b_e+2,a_s-2:a_e+2)
    integer, intent(in)    :: side
    integer, intent(in)    :: direction

    integer :: a
    integer :: a_s_buf, a_e_buf
    integer :: b_s_buf, b_e_buf
    integer :: c_s_buf, c_e_buf
    integer :: na, nb, nc, nn

    call get_buffer_indices(side, direction, .true., &
                            a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)

    ibuf = ibuf + 1
    nn = nb * nc * n1 * n2 * n3 * n4

    !$omp parallel do
    do a = a_s_buf, a_e_buf
      x(1:n4,1:n3,1:n2,1:n1,c_s_buf:c_e_buf,b_s_buf:b_e_buf,a) = &
          RESHAPE(rbuf(ibuf+(a-a_s_buf)*nn:ibuf+(a-a_s_buf+1)*nn-1), (/n4,n3,n2,n1,nc,nb/))
    enddo
    !$omp end parallel do

    ibuf = ibuf + na*nn - 1

  end subroutine read_buffer_4d

  subroutine get_buffer_indices(side, direction, ghost, a_s_buf, a_e_buf, b_s_buf, b_e_buf, c_s_buf, c_e_buf, na, nb, nc)
    integer, intent(in)  :: side
    integer, intent(in)  :: direction
    logical, intent(in)  :: ghost
    integer, intent(out) :: a_s_buf
    integer, intent(out) :: a_e_buf
    integer, intent(out) :: b_s_buf
    integer, intent(out) :: b_e_buf
    integer, intent(out) :: c_s_buf
    integer, intent(out) :: c_e_buf
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nc

    if (side == 0) then
      if (direction == 1) then
        if (ghost) then
          a_s_buf = a_s - 2
          a_e_buf = a_s - 1
        else
          a_s_buf = a_s
          a_e_buf = a_s + 1
        endif
        b_s_buf = b_s - 2
        b_e_buf = b_e + 2
        c_s_buf = c_s - 2
        c_e_buf = c_e + 2
      else if (direction == 2) then
        a_s_buf = a_s - 2
        a_e_buf = a_e + 2
        if (ghost) then
          b_s_buf = b_s - 2
          b_e_buf = b_s - 1
        else
          b_s_buf = b_s
          b_e_buf = b_s + 1
        endif
        c_s_buf = c_s - 2
        c_e_buf = c_e + 2
      else if (direction == 3) then
        a_s_buf = a_s - 2
        a_e_buf = a_e + 2
        b_s_buf = b_s - 2
        b_e_buf = b_e + 2
        if (ghost) then
          c_s_buf = c_s - 2
          c_e_buf = c_s - 1
        else
          c_s_buf = c_s
          c_e_buf = c_s + 1
        endif
      endif
    else if (side == 1) then
      if (direction == 1) then
        if (ghost) then
          a_s_buf = a_e + 1
          a_e_buf = a_e + 2
        else
          a_s_buf = a_e - 1
          a_e_buf = a_e
        endif
        b_s_buf = b_s - 2
        b_e_buf = b_e + 2
        c_s_buf = c_s - 2
        c_e_buf = c_e + 2
      else if (direction == 2) then
        a_s_buf = a_s - 2
        a_e_buf = a_e + 2
        if (ghost) then
          b_s_buf = b_e + 1
          b_e_buf = b_e + 2
        else
          b_s_buf = b_e - 1
          b_e_buf = b_e
        endif
        c_s_buf = c_s - 2
        c_e_buf = c_e + 2
      else if (direction == 3) then
        a_s_buf = a_s - 2
        a_e_buf = a_e + 2
        b_s_buf = b_s - 2
        b_e_buf = b_e + 2
        if (ghost) then
          c_s_buf = c_e + 1
          c_e_buf = c_e + 2
        else
          c_s_buf = c_e - 1
          c_e_buf = c_e
        endif
      endif
    endif

    na = a_e_buf - a_s_buf + 1
    nb = b_e_buf - b_s_buf + 1
    nc = c_e_buf - c_s_buf + 1

  end subroutine get_buffer_indices

  subroutine sum_mpi(x)
    real, intent(inout) :: x

#ifdef MPI_BT
    real    :: xsum
    integer :: ierr

    call MPI_Allreduce(x, xsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    x = xsum
#else
    ! Prevent warning if MPI is not enabled
    if (x > 0) continue
#endif

  end subroutine sum_mpi

  subroutine bcast_mpi(x)
    real, intent(inout) :: x

#ifdef MPI_BT
    integer :: ierr
    call MPI_Bcast(x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
    ! Prevent warning if MPI is not enabled
    if (x > 0) continue
#endif

  end subroutine bcast_mpi

  subroutine send_mpi(x, dst, tag)
    real,     intent(in) :: x
    integer,  intent(in) :: dst
    integer,  intent(in) :: tag

#ifdef MPI_BT
    integer :: ierr
    call MPI_Send(x, 1, MPI_DOUBLE_PRECISION, dst, tag, MPI_COMM_WORLD, ierr)
#else
    if (x   > 0) continue
    if (dst > 0) continue
    if (tag > 0) continue
#endif

  end subroutine send_mpi

  subroutine recv_mpi(x, tag)
    real,     intent(out) :: x
    integer,  intent(in)  :: tag

#ifdef MPI_BT
integer :: mpistat(MPI_STATUS_SIZE)
    integer :: ierr
    call MPI_Recv(x, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, mpistat, ierr)
#else
    x = 0.0
    if (tag > 0) continue
#endif

  end subroutine recv_mpi

end module rice_mpi
