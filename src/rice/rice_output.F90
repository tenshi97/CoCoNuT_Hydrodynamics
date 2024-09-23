module rice_output
  use hdf5
  use rice_config,  only: nr, ntheta, nphi, nmu, npsi, neps, nflav
  use rice_grid,    only: a_s, a_e, b_s, b_e, c_s, c_e
  implicit none

  private

  public :: write_output, read_output, init_hdf5, finish_hdf5

#ifdef MPI_BT
  logical, parameter :: compression = .false.
#else
  logical, parameter :: compression = .true.
#endif
  integer, parameter :: compression_level = 9

  interface write_output
   module procedure write_output_int, write_output_filename
  end interface

  interface read_output
   module procedure read_output_int, read_output_filename
  end interface

  interface file_read
   module procedure read_scalar, &
                    read_scalar_int, &
                    read_array_3d, &
                    read_array_4d, &
                    read_array_7d
  end interface

  interface file_write
   module procedure write_scalar, &
                    write_scalar_int, &
                    write_array_1d, &
                    write_array_2d, &
                    write_array_3d, &
                    write_array_4d, &
                    write_array_5d, &
                    ! write_array_6d, &
                    write_array_7d
  end interface

contains

  subroutine check_error(error)
    integer, intent(in) :: error

    if (error < 0) stop 'HDF5 error'

  end subroutine check_error

  subroutine init_hdf5
    integer :: error
    call h5open_f(error)
    call check_error(error)
  end subroutine init_hdf5

  subroutine finish_hdf5
    integer :: error
    call h5close_f(error)
    call check_error(error)
  end subroutine finish_hdf5

  subroutine read_output_filename(filename)
    use rice_mpi
    use rice_config, only: ic_reset
    use rice_grid,   only: f, vfluid_old, alpha_old, beta_old, phiconf_old, time, index, nstep

    character(*), intent(in) :: filename

    integer(HID_T)    :: file_id, prop_id
    integer           :: error

    integer :: ext_7d(7*2), ext_7d_l(7*2), ext_7d_g(7*2)

    ext_7d   = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, 1,     nphi,  1,     ntheta, 1,     nr   /)
    ext_7d_l = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, c_s,   c_e,   b_s,   b_e,    a_s,   a_e  /)
    ext_7d_g = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, c_s-2, c_e+2, b_s-2, b_e+2,  a_s-2, a_e+2/)

    print*, '(output) reading ', trim(filename)

    ! File access property list
    call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! MPI communicator information
    call h5pset_fapl_mpio_f(prop_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    call check_error(error)
#endif

    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp=prop_id)
    call check_error(error)

    call h5pclose_f(prop_id, error)
    call check_error(error)

    call file_read(index,       'index',   file_id)
    call file_read(nstep,       'nstep',   file_id)
    call file_read(time,        'time',    file_id)

    call file_read(f,'f', file_id, ext_7d_g, ext_7d_l, ext_7d)
    call file_read(vfluid_old  (1:3,                           1:nphi,1:ntheta,1:nr), 'vfluid',  file_id)
    call file_read(alpha_old   (                               1:nphi,1:ntheta,1:nr), 'alpha',   file_id)
    call file_read(beta_old    (1:3,                           1:nphi,1:ntheta,1:nr), 'beta',    file_id)
    call file_read(phiconf_old (                               1:nphi,1:ntheta,1:nr), 'phiconf', file_id)

    if (ic_reset) call reset_counters

    call h5fclose_f(file_id, error)
    call check_error(error)

  end subroutine read_output_filename

  subroutine read_output_int(index)
    integer,  intent(in) :: index

    character(len=200) :: filename

    write(filename, "('output/snap_', I6.6, '.h5')") index

    call read_output_filename(filename)

  end subroutine read_output_int

  subroutine read_array_3d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 3
    real,           intent(out) :: x(:,:,:)
    character(*),   intent(in)  :: name
    integer(HID_T), intent(in)  :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed out
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be read to
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call read_array_nd(name, ndims, file_id, &
    x3=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine read_array_3d

  subroutine read_array_4d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 4
    real,           intent(out) :: x(:,:,:,:)
    character(*),   intent(in)  :: name
    integer(HID_T), intent(in)  :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed out
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be read to
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call read_array_nd(name, ndims, file_id, &
    x4=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine read_array_4d

  subroutine read_array_7d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 7
    real,           intent(out) :: x(:,:,:,:,:,:,:)
    character(*),   intent(in)  :: name
    integer(HID_T), intent(in)  :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed out
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be read to
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call read_array_nd(name, ndims, file_id, &
    x7=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine read_array_7d

  subroutine read_array_nd(name, ndims, file_id, &
                           x1, x2, x3, x4, x5, x6, x7, &
                           ext_array, ext_local, ext_global)

    character(*),      intent(in)  :: name
    integer,           intent(in)  :: ndims
    integer(HID_T),    intent(in)  :: file_id
    real,    optional, intent(out) :: x1(:)
    real,    optional, intent(out) :: x2(:,:)
    real,    optional, intent(out) :: x3(:,:,:)
    real,    optional, intent(out) :: x4(:,:,:,:)
    real,    optional, intent(out) :: x5(:,:,:,:,:)
    real,    optional, intent(out) :: x6(:,:,:,:,:,:)
    real,    optional, intent(out) :: x7(:,:,:,:,:,:,:)
    integer, optional, intent(in)  :: ext_array (ndims*2) ! extents of the array being passed out
    integer, optional, intent(in)  :: ext_local (ndims*2) ! subset of the array to be read to
    integer, optional, intent(in)  :: ext_global(ndims*2) ! extents of the full grid

    integer(HSIZE_T)    :: shape_array(ndims), shape_local(ndims), shape_global(ndims), shapex(ndims)
    integer(HSIZE_T)    :: offset_array(ndims), offset_global(ndims)
    integer(HID_T)      :: dset_id, dspace_id, mspace_id, prop_id
    integer             :: error
    logical             :: is_simple
    integer             :: ndims_data
    integer(HSIZE_T)    :: shape_data(ndims), maxshape_data(ndims)

    if (present(ext_array)) then
      shape_array  = ext_array (2::2) - ext_array (1::2) + 1 ! shape of array being passed in
      shape_local  = ext_local (2::2) - ext_local (1::2) + 1 ! shape of array to write
      shape_global = ext_global(2::2) - ext_global(1::2) + 1 ! shape of the full global array

      offset_array  = ext_local(1::2) - ext_array (1::2) ! offset of array to write against the array being passed in
      offset_global = ext_local(1::2) - ext_global(1::2) ! offset of array to write against the full global array
    else
      if (present(x1)) then
        shapex = shape(x1)
      else if (present(x2)) then
        shapex = shape(x2)
      else if (present(x3)) then
        shapex = shape(x3)
      else if (present(x4)) then
        shapex = shape(x4)
      else if (present(x5)) then
        shapex = shape(x5)
      else if (present(x6)) then
        shapex = shape(x6)
      else if (present(x7)) then
        shapex = shape(x7)
      endif
      shape_array  = shapex
      shape_local  = shapex
      shape_global = shapex

      offset_array (:) = 0
      offset_global(:) = 0
    endif

    ! Open dataset
    call h5dopen_f(file_id, name, dset_id, error)
    call check_error(error)

    ! Get dataspace
    call h5dget_space_f(dset_id, dspace_id, error)
    call check_error(error)

    ! Check that the space is simple
    call h5sis_simple_f(dspace_id, is_simple, error)
    call check_error(error)
    if (.not. is_simple) stop "Dataspace is not simple"

    ! Check number of dimensions
    call h5sget_simple_extent_ndims_f(dspace_id, ndims_data, error)
    call check_error(error)
    if (ndims_data /= ndims) stop "Dimensions of file data do not match"

    ! Check shape
    call h5sget_simple_extent_dims_f(dspace_id, shape_data, maxshape_data, error)
    call check_error(error)
    if (.not. all(   shape_data == shape_global)) stop "Shape of file data does not match"
    if (.not. all(maxshape_data == shape_global)) stop "Shape of file data does not match"

    ! Use hyperslab to select the part of the file to read from
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_global, shape_local, error)
    call check_error(error)

    ! Define dataset in memory
    call h5screate_simple_f(ndims, shape_array, mspace_id, error)
    call check_error(error)

    ! Use hyperslab to select the part of the array to write to
    call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, offset_array, shape_local, error)
    call check_error(error)

    ! Create property list for collective dataset read
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective read
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Read file
    if (present(x1)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x1, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x2)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x2, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x3)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x3, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x4)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x4, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x5)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x5, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x6)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x6, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x7)) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x7, shape_local, error, mspace_id, dspace_id, prop_id)
    endif
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

    ! Close memory space
    call h5sclose_f(mspace_id, error)
    call check_error(error)

    ! Close dataspace
    call h5sclose_f(dspace_id, error)
    call check_error(error)

    ! Close property list
    call h5pclose_f(prop_id, error)
    call check_error(error)

  end subroutine read_array_nd

  subroutine read_scalar(x, name, file_id)
    real,           intent(out) :: x
    character(*),   intent(in)  :: name
    integer(HID_T), intent(in)  :: file_id

    integer(HSIZE_T)    :: xshape(1) = 0
    integer(HID_T)      :: dset_id, prop_id
    integer             :: error

    ! Open dataset
    call h5dopen_f(file_id, name, dset_id, error)
    call check_error(error)

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective write
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Read file
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error, xfer_prp=prop_id)
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

  end subroutine read_scalar

  subroutine read_scalar_int(x, name, file_id)
    integer,        intent(out) :: x
    character(*),   intent(in)  :: name
    integer(HID_T), intent(in)  :: file_id

    integer(HSIZE_T)    :: xshape(1) = 0
    integer(HID_T)      :: dset_id, prop_id
    integer             :: error

    ! Open dataset
    call h5dopen_f(file_id, name, dset_id, error)
    call check_error(error)

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective write
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Read file
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error, xfer_prp=prop_id)
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

  end subroutine read_scalar_int

  subroutine write_output_filename(filename)
    use rice_mpi
    character(200), intent(in) :: filename

    integer(HID_T) :: file_id, prop_id
    integer :: error

    print*, '(output) writing ', trim(filename)

    call execute_command_line('mkdir -p output/')

    ! File access property list
    call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! MPI communicator information
    call h5pset_fapl_mpio_f(prop_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    call check_error(error)
#endif

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp=prop_id)
    call check_error(error)

    call h5pclose_f(prop_id, error)
    call check_error(error)

    call write_file(file_id)

    call h5fclose_f(file_id, error)
    call check_error(error)

  end subroutine write_output_filename

  subroutine write_output_int
    use rice_grid, only: index
    character(len=200) :: filename

    index = index + 1
    write(filename, "('output/snap_', I6.6, '.h5')") index

    call write_output_filename(filename)

  end subroutine write_output_int

  subroutine write_file(file_id)
    use rice_config, only: clight, ic_reset, cartoon_grid
    use rice_grid, only: index, nstep, f, mu, psi, r, r_if, theta, theta_if, phi, phi_if, &
                         domega, f_eq, kappa_a, kappa_s, eps, eps_if, vfluid, &
                         area_r, area_theta, area_phi, volume, time, &
                         area_r_spherical, area_theta_spherical, area_phi_spherical, &
                         volume_spherical, &
                         alpha, beta, phiconf, u_com

    ! extents: global       local          local with ghosts
    integer :: ext_3d(3*2), ext_3d_l(3*2), ext_3d_g(3*2)
    integer :: ext_4d(4*2), ext_4d_l(4*2), ext_4d_g(4*2)
    integer :: ext_5d(5*2), ext_5d_l(5*2), ext_5d_g(5*2)
    integer :: ext_7d(7*2), ext_7d_l(7*2), ext_7d_g(7*2)

    integer(HID_T), intent(in) :: file_id

    ext_3d   = (/                                       1,     nphi,  1,     ntheta, 1,     nr   /)
    ext_3d_l = (/                                       c_s,   c_e,   b_s,   b_e,    a_s,   a_e  /)
    ext_3d_g = (/                                       c_s-2, c_e+2, b_s-2, b_e+2,  a_s-2, a_e+2/)

    ext_4d   = (/1, 3,                                  1,     nphi,  1,     ntheta, 1,     nr   /)
    ext_4d_l = (/1, 3,                                  c_s,   c_e,   b_s,   b_e,    a_s,   a_e  /)
    ext_4d_g = (/1, 3,                                  c_s-2, c_e+2, b_s-2, b_e+2,  a_s-2, a_e+2/)

    ext_5d   = (/1, nflav, 1, neps,                     1,     nphi,  1,     ntheta, 1,     nr   /)
    ext_5d_l = (/1, nflav, 1, neps,                     c_s,   c_e,   b_s,   b_e,    a_s,   a_e  /)
    ext_5d_g = (/1, nflav, 1, neps,                     c_s-2, c_e+2, b_s-2, b_e+2,  a_s-2, a_e+2/)

    ext_7d   = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, 1,     nphi,  1,     ntheta, 1,     nr   /)
    ext_7d_l = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, c_s,   c_e,   b_s,   b_e,    a_s,   a_e  /)
    ext_7d_g = (/1, nflav, 1, neps, 1, npsi, -nmu, nmu, c_s-2, c_e+2, b_s-2, b_e+2,  a_s-2, a_e+2/)

    if (ic_reset) then
      call file_write(0,          'index',      file_id)
      call file_write(0,          'nstep',      file_id)
      call file_write(0.0,        'time',       file_id)
    else
      call file_write(index,      'index',      file_id)
      call file_write(nstep,      'nstep',      file_id)
      call file_write(time,       'time',       file_id)
    endif

    ! Only write real cells, not ghost cells
    call file_write(f,        'f',        file_id, ext_7d_g,         ext_7d_l,      ext_7d       )
    call file_write(f_eq,     'f_eq',     file_id, ext_5d_g,         ext_5d_l,      ext_5d       )
    call file_write(kappa_a,  'kappa_a',  file_id, ext_5d_g,         ext_5d_l,      ext_5d       )
    call file_write(kappa_s,  'kappa_s',  file_id, ext_5d_g,         ext_5d_l,      ext_5d       )
    call file_write(r,        'r',        file_id, (/-1, nr+2/),     (/1, nr/),     (/1, nr/)    )
    call file_write(r_if,     'r_if',     file_id, (/-1, nr+1/),     (/0, nr/),     (/0, nr/)    )
    call file_write(theta,    'theta',    file_id, (/-1, ntheta+2/), (/1, ntheta/), (/1, ntheta/))
    call file_write(theta_if, 'theta_if', file_id, (/-1, ntheta+1/), (/0, ntheta/), (/0, ntheta/))
    call file_write(phi,      'phi',      file_id, (/-1, nphi+2/),   (/1, nphi/),   (/1, nphi/)  )
    call file_write(phi_if,   'phi_if',   file_id, (/-1, nphi+1/),   (/0, nphi/),   (/0, nphi/)  )
    call file_write(vfluid,   'vfluid',   file_id, ext_4d_g,         ext_4d_l,      ext_4d       )
    call file_write(alpha,    'alpha',    file_id, ext_3d_g,         ext_3d_l,      ext_3d       )
    call file_write(beta,     'beta',     file_id, ext_4d_g,         ext_4d_l,      ext_4d       )
    call file_write(phiconf,  'phiconf',  file_id, ext_3d_g,         ext_3d_l,      ext_3d       )

    call file_write(u_com,  'u_com',  file_id)
    call file_write(mu,     'mu',     file_id)
    call file_write(psi,    'psi',    file_id)
    call file_write(eps,    'eps',    file_id)
    call file_write(eps_if, 'eps_if', file_id)
    call file_write(domega, 'domega', file_id)
    call file_write(clight, 'clight', file_id)


    if (cartoon_grid) then
      call file_write(area_r_spherical,                 'area_r',               file_id)
      call file_write(area_theta_spherical,             'area_theta',           file_id)
      call file_write(area_phi_spherical,               'area_phi',             file_id)
      call file_write(area_r,                           'area_r_cartoon',       file_id)
      call file_write(area_theta,                       'area_theta_cartoon',   file_id)
      call file_write(area_phi,                         'area_phi_cartoon',     file_id)
      call file_write(volume_spherical(1:nphi,1:ntheta,1:nr), 'volume',         file_id)
      call file_write(volume          (1:nphi,1:ntheta,1:nr), 'volume_cartoon', file_id)
    else
      call file_write(area_r,                           'area_r',             file_id)
      call file_write(area_theta,                       'area_theta',         file_id)
      call file_write(area_phi,                         'area_phi',           file_id)
      call file_write(volume(1:nphi,1:ntheta,1:nr),     'volume',             file_id)
    endif

  end subroutine write_file

  subroutine write_scalar(x, name, file_id)
    real,           intent(in) :: x
    character(*),   intent(in) :: name
    integer(HID_T), intent(in) :: file_id

    integer(HSIZE_T), parameter  :: xshape(1) = 0
    integer(HID_T)    :: dspace_id, dset_id, prop_id
    integer           :: error

    ! Create dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    call check_error(error)

    ! Create dataset in file
    call h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    call check_error(error)

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective write
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Write to file
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error, xfer_prp=prop_id)
    call check_error(error)

    ! Close property list
    call h5pclose_f(prop_id, error)
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

    ! Closet dataspace
    call h5sclose_f(dspace_id, error)
    call check_error(error)

  end subroutine write_scalar

  subroutine write_scalar_int(x, name, file_id)
    integer,        intent(in) :: x
    character(*),   intent(in) :: name
    integer(HID_T), intent(in) :: file_id

    integer(HSIZE_T), parameter  :: xshape(1) = 0
    integer(HID_T)    :: dspace_id, dset_id, prop_id
    integer           :: error

    ! Create dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    call check_error(error)

    ! Create dataset in file
    call h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    call check_error(error)

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective write
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Write to file
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error, xfer_prp=prop_id)
    call check_error(error)

    ! Close property list
    call h5pclose_f(prop_id, error)
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

    ! Closet dataspace
    call h5sclose_f(dspace_id, error)
    call check_error(error)

  end subroutine write_scalar_int

  subroutine write_array_1d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 1
    real,              intent(in) :: x(:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x1=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_1d

  subroutine write_array_2d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 2
    real,              intent(in) :: x(:,:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x2=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_2d

  subroutine write_array_3d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 3
    real,              intent(in) :: x(:,:,:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x3=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_3d

  subroutine write_array_4d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 4
    real,              intent(in) :: x(:,:,:,:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x4=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_4d

  subroutine write_array_5d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 5
    real,              intent(in) :: x(:,:,:,:,:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x5=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_5d

  ! subroutine write_array_6d(x, name, file_id, ext_array, ext_local, ext_global)
  !   integer, parameter  :: ndims = 6
  !   real,              intent(in) :: x(:,:,:,:,:,:)
  !   character(*),      intent(in) :: name
  !   integer(HID_T),    intent(in) :: file_id
  !   integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
  !   integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
  !   integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

  !   call write_array_nd(name, ndims, file_id, &
  !   x6=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  ! end subroutine write_array_6d

  subroutine write_array_7d(x, name, file_id, ext_array, ext_local, ext_global)
    integer, parameter  :: ndims = 7
    real,              intent(in) :: x(:,:,:,:,:,:,:)
    character(*),      intent(in) :: name
    integer(HID_T),    intent(in) :: file_id
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    call write_array_nd(name, ndims, file_id, &
    x7=x, ext_array=ext_array, ext_local=ext_local, ext_global=ext_global)

  end subroutine write_array_7d

  subroutine write_array_nd(name, ndims, file_id, &
                            x1, x2, x3, x4, x5, x6, x7, &
                            ext_array, ext_local, ext_global)
    use rice_mpi, only: id
    character(*),      intent(in) :: name
    integer,           intent(in) :: ndims
    integer(HID_T),    intent(in) :: file_id
    real,    optional, intent(in) :: x1(:)
    real,    optional, intent(in) :: x2(:,:)
    real,    optional, intent(in) :: x3(:,:,:)
    real,    optional, intent(in) :: x4(:,:,:,:)
    real,    optional, intent(in) :: x5(:,:,:,:,:)
    real,    optional, intent(in) :: x6(:,:,:,:,:,:)
    real,    optional, intent(in) :: x7(:,:,:,:,:,:,:)
    integer, optional, intent(in) :: ext_array (ndims*2) ! extents of the array being passed in
    integer, optional, intent(in) :: ext_local (ndims*2) ! subset of the array to be written
    integer, optional, intent(in) :: ext_global(ndims*2) ! extents of the full grid

    integer(HSIZE_T)    :: shape_array(ndims), shape_local(ndims), shape_global(ndims), shapex(ndims)
    integer(HSIZE_T)    :: offset_array(ndims), offset_global(ndims)
    integer(HSIZE_T)    :: chunk(ndims)
    integer(HID_T)      :: dspace_id, dset_id, prop_id, mspace_id
    integer             :: error

    logical :: parallel_write

    if (present(ext_array)) then
      shape_array  = ext_array (2::2) - ext_array (1::2) + 1 ! shape of array being passed in
      shape_local  = ext_local (2::2) - ext_local (1::2) + 1 ! shape of array to write
      shape_global = ext_global(2::2) - ext_global(1::2) + 1 ! shape of the full global array

      offset_array  = ext_local(1::2) - ext_array (1::2) ! offset of array to write against the array being passed in
      offset_global = ext_local(1::2) - ext_global(1::2) ! offset of array to write against the full global array
    else
      if (present(x1)) then
        shapex = shape(x1)
      else if (present(x2)) then
        shapex = shape(x2)
      else if (present(x3)) then
        shapex = shape(x3)
      else if (present(x4)) then
        shapex = shape(x4)
      else if (present(x5)) then
        shapex = shape(x5)
      else if (present(x6)) then
        shapex = shape(x6)
      else if (present(x7)) then
        shapex = shape(x7)
      endif
      shape_array  = shapex
      shape_local  = shapex
      shape_global = shapex

      offset_array (:) = 0
      offset_global(:) = 0
    endif

    if (all(shape_local == shape_global)) then
      parallel_write = .false.
    else
      parallel_write = .true.
    endif

    ! Create dataspace
    call h5screate_simple_f(ndims, shape_global, dspace_id, error)
    call check_error(error)

    ! Create the dataset creation property list, add the gzip
    ! compression filter and set the chunk size.
    call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, error)
    call check_error(error)
    if (compression) then
      chunk = shape_local
      call h5pset_deflate_f(prop_id, compression_level, error)
      call check_error(error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check_error(error)
    endif

    ! Create dataset in file
    call h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error, prop_id)
    call check_error(error)

    ! Close property list
    call h5pclose_f(prop_id, error)
    call check_error(error)

    ! Use hyperslab to select the part of the file to write to
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset_global, shape_local, error)
    call check_error(error)

    ! If not a parallel write, then only master writes
    if ((.not. parallel_write) .and. (id /= 0)) then
      call h5sselect_none_f(dspace_id, error)
      call check_error(error)
    endif

    ! Define dataset in memory
    call h5screate_simple_f(ndims, shape_array, mspace_id, error)
    call check_error(error)

    ! Use hyperslab to select the part of the array to write
    call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, offset_array, shape_local, error)
    call check_error(error)

    ! Select none if not writing
    if ((.not. parallel_write) .and. (id /= 0)) then
      call h5sselect_none_f(mspace_id, error)
      call check_error(error)
    endif

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, error)
    call check_error(error)

#ifdef MPI_BT
    ! Set collective write
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error(error)
#endif

    ! Write to file
    if (present(x1)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x1, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x2)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x2, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x3)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x3, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x4)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x4, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x5)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x5, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x6)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x6, shape_local, error, mspace_id, dspace_id, prop_id)
    else if (present(x7)) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x7, shape_local, error, mspace_id, dspace_id, prop_id)
    endif
    call check_error(error)

    ! Close dataset
    call h5dclose_f(dset_id, error)
    call check_error(error)

    ! Close dataspace
    call h5sclose_f(dspace_id, error)
    call check_error(error)

    ! Close memory space
    call h5sclose_f(mspace_id, error)
    call check_error(error)

    ! Close property list
    call h5pclose_f(prop_id, error)
    call check_error(error)

  end subroutine write_array_nd

  subroutine reset_counters
    use rice_grid, only: index, nstep, time

    index = 0
    nstep = 0
    time  = 0.0
  end subroutine reset_counters

end module rice_output
