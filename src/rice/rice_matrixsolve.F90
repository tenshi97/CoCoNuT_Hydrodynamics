#define ANALYTIC_SCATTERING

module rice_matrixsolve
  implicit none

contains

#ifndef ANALYTIC_SCATTERING
  subroutine solve_matrix(matrix_size, matrix, rhs)
    integer, intent(in)     :: matrix_size
    real,    intent(inout)  :: matrix(matrix_size,matrix_size)
    real,    intent(inout)  :: rhs(matrix_size)

    integer, parameter :: nrhs=1
    integer :: lda, ldb
    integer :: info

    integer :: ipiv(matrix_size)

    lda = matrix_size
    ldb = matrix_size

    call dgesv(matrix_size, nrhs, matrix, lda, ipiv, rhs, ldb, info)

    if (info > 0) then
      print*,'A is singular, solution could not be computed'
      stop
    endif

  end subroutine solve_matrix
#endif

end module rice_matrixsolve
