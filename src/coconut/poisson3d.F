#if (DIMN==3) && !defined(ENFORCE_SPHERICAL_METRIC)
c     ================================================================

      MODULE poisson3d

c     ================================================================
c
c     AUTHOR: Bernhard Mueller
c     DATE:   3 April 2013
c
c     This module solves the Poisson equation in 3D spherical
c     polar coordinates using a decomposition into eigenfunctions of the
c     angular component of the discrete Laplacian. The transformation
c     between the real and spectral domains is handled by FFTs for the
c     phi-direction and MVMs for the theta-direction.
c
c     Input:
c     nr:   Number of radial zones
c     nth:  Number of zones in theta-direction
c     nphi: Number of zones in r-direction
c     rho:  Source term in the Poisson equation
c     
c     Output:
c     phi:  Solution field
c
c     ----------------------------------------------------------------

      USE precision

      IMPLICIT NONE

      SAVE

c     ----------------------------------------------------------------
c     Shared module variables
c     ----------------------------------------------------------------
c     For both the FFT and the poisson solver:
      INTEGER (KIND=ik) :: nr, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc

c     ----------------------------------------------------------------
c     For the FFT
      INTEGER (KIND=ik) :: logn
      REAL (KIND=rk) :: pi

      INTEGER (KIND=ik) :: binary (0:20)
      INTEGER (KIND=ik), ALLOCATABLE :: index (:)


c     For the Poisson solver

      INTEGER :: nth_proc

      REAL (KIND=rk) :: dtheta, dphi
      REAL (KIND=rk) :: lambda_phi

      REAL (KIND=rk), ALLOCATABLE :: r_if (:)

      REAL (KIND=rk), ALLOCATABLE :: a (:), b(:), c (:)
      REAL (KIND=rk), ALLOCATABLE :: lap_th (:,:)
      REAL (KIND=rk), ALLOCATABLE :: vol_th (:,:)
      
      REAL (KIND=rk), ALLOCATABLE :: matrix (:,:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: lambda (:,:)

      REAL (KIND=rk), ALLOCATABLE :: da_r (:), dv_r (:)
      REAL (KIND=rk), ALLOCATABLE :: diag0 (:), offdiag0 (:)

      CONTAINS


c     ================================================================

      SUBROUTINE setup_poisson_and_fft (m_in, n_in, o_in,
     &     n_s_in, n_e_in, o_e_in, o_s_in, r_if_in, theta_if, theta)

c     ================================================================

      USE precision

      IMPLICIT NONE

c     ----------------------------------------------------------------
c     Subroutine arguments:
      INTEGER (KIND=ik), INTENT (IN) :: m_in, n_in, o_in,
     &     n_s_in, n_e_in, o_e_in, o_s_in
      REAL (KIND=rk), INTENT (IN) :: r_if_in (0:m_in)
      REAL (KIND=rk), INTENT (IN) :: theta_if (0:n_in)
      REAL (KIND=rk), INTENT (IN) :: theta    (1:n_in)

c     ----------------------------------------------------------------
      REAL (KIND=rk) :: r (1:m_in)
      REAL (KIND=rk) :: w (n_in/2)
      REAL (KIND=rk) :: vec (n_in/2,n_in/2)
      REAL (KIND=rk) :: work (3*n_in) 

      INTEGER (KIND=ik) :: i, j, jj, k, kk, istat, nth, nth2
      INTEGER (KIND=ik) :: info
c     ----------------------------------------------------------------

      pi = ACOS (-1.0d0)
      PRINT*,'Dimensions',m_in,n_in,o_in,n_s_in,n_e_in,o_s_in,o_e_in
      PRINT*,r_if_in(:)
      PRINT*,theta_if(:)
      nr  = m_in
      n   = n_in
      o   = o_in
      n_s = n_s_in
      n_e = n_e_in
      o_s = o_s_in
      o_e = o_e_in

      n_loc = n_e - n_s + 1
      o_loc = o_e - o_s + 1

c     -----------------------------------------------------------------
c     FFT setup:
c     Determine log_2 (qz), integer arithmetic
      j = o
      logn = 0
      DO i = 1, o
         j = j / 2
         IF (j .EQ. 0) EXIT
         logn = logn + 1
         PRINT *,'log_2 (qz) =',logn
      END DO

      ALLOCATE (index(0:o-1), STAT=istat)

      DO i = 0, o-1
         j = i
         DO k = logn-1, 0, -1
            binary (k) = j / 2 ** k
            j = j - binary (k) * 2 ** k
         END DO

         j=0
         DO k = 0, logn-1
            j = j + binary (logn-1-k) * 2 ** k
         END DO
         index (i) = j
      END DO

c     -----------------------------------------------------------------
c     Setup Poisson solver
      nth = n
      nth2 = n / 2
      nth_proc = n / n_loc

      ALLOCATE (lap_th (2,1:nth/2), STAT = istat)
      ALLOCATE (vol_th (1,1:nth), STAT = istat)
      ALLOCATE (a (1:nth), b (1:nth), c (1:nth), STAT = istat)
      ALLOCATE (matrix (1:nth2,1:nth2,o_s:o_e,0:1) ,STAT = istat)
      ALLOCATE (lambda (1:nth,o_s:o_e), STAT = istat)

      ALLOCATE (diag0 (1:nr), offdiag0 (1:nr-1), STAT = istat)
      ALLOCATE (da_r (0:nr), dv_r (1:nr), STAT = istat)
      ALLOCATE (r_if (0:nr), STAT = istat)

      r_if (0:nr) = r_if_in (0:nr)

      DO i = 0, nr
         da_r (i) = r_if (i) ** 2
      END DO
      DO i = 1, nr
         r (i) = 0.5d0 * (r_if (i) + r_if (i-1))
         dv_r (i) = (r_if (i) ** 3 - r_if (i-1) ** 3) / 3.0d0
      END DO      
      DO i = 1, nr-1
         offdiag0 (i) = - da_r (i) / (r (i+1) - r(i))
         diag0 (i) = - offdiag0 (i)
      END DO
      diag0 (nr) = 0.0d0 
      DO i = 2, nr
         diag0 (i) = diag0 (i) - offdiag0 (i-1)
      END DO

c     equidistant grid in theta assumed at the moment
      dtheta = pi / REAL (n, KIND=rk)
      dphi = 2.0d0 * pi / REAL (o, KIND=rk)
      

c     -----------------------------------------------------------------
c     Compute eigenfunctions
c     -----------------------------------------------------------------

      DO kk = o_s, o_e
         k = index (kk-1)         
c     k runs over frequency indices stored by local task
         
c         lambda_phi = 
c     &        (2.0d0 *  SIN (0.5d0 * DBLE (k) * dphi) / dphi) ** 2 *
c     &        (0.5d0 * dtheta / SIN ( 0.5d0 * dtheta))
         lambda_phi = 
     &        (2.0d0 *  SIN (0.5d0 * DBLE (k) * dphi) / dphi) ** 2 *
c     &        (0.5d0 * dtheta / SIN ( 0.5d0 * dtheta))
     &        dtheta         


         DO j = 1, nth
            a (j) = SIN (theta_if (j-1)) / dtheta ** 2
            b (j) = - (SIN (theta_if (j-1)) + SIN (theta_if (j))) /
c     &           dtheta ** 2 - lambda_phi / SIN (theta (j))
     &           dtheta ** 2 - lambda_phi / SIN (theta (j)) / dtheta 
c     &           * (COS (theta_if (j-1)) - COS (theta_if (j))) / dtheta
            c (j) = SIN (theta_if (j)) / dtheta ** 2
c            a (j) = SIN (theta_if (j-1)) / dtheta ** 1
c            b (j) = - (SIN (theta_if (j-1)) + SIN (theta_if (j))) /
c     &           dtheta ** 1 - lambda_phi / SIN (theta (j))
c            c (j) = SIN (theta_if (j)) / dtheta ** 1
         END DO

         lap_th (:,:) = 0.0d0
         vol_th (:,:) = 0.0d0

c     Eigenvectors with EVEN parity
         DO j = 1, nth2
            lap_th (1,j) = b (j)
            IF (j .LE. nth2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth2) lap_th (1,j) = b (j) + c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
c            vol_th (1,j) = 
c     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
c     &           dtheta**0
            PRINT *,'Poisson-3D (operator matrix)',j,lap_th(:,j),
     &           vol_th(1,j)
         END DO
         
         CALL DSBGV ('V','L',nth2,1,0,lap_th(:,1:nth2),2,
     &        vol_th(1,1:nth2),1,w,vec,nth2,work,info)
         
         DO j = 1, nth
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
         END DO
         DO j = 1, nth2 
            matrix (1:nth2,j,kk,0) = vec(1:nth2,j)
         END DO
         lambda (1:nth2,kk) = DBLE (w (1:nth2))
         IF (info .NE. 0) STOP 'DSBGV failed!'
         
c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth/2
            DO i = 1, nth/2
               PRINT *, k, vec (i,j), 
     &              j, i, COS (theta (i)),'eigenvec (even)'
            END DO
         END DO

c     effective l*(l+1)
         DO j = 1, nth/2
            PRINT *,'l*(l+1), even',j, k, w (j), lambda_phi
         END DO

c     Eigenvectors with ODD parity
         lap_th (:,:) = 0.0d0

         DO j = 1, nth2
            lap_th (1,j) = b (j)
            IF (j .LE. nth2-1) lap_th (2,j) = c (j)
            IF (j .EQ. nth2) lap_th (1,j) = b (j) - c (j)
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
c            vol_th (1,j) = 
c     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
c     &           dtheta**0
            PRINT *,'Poisson-3D (operator matrix)',j,lap_th(:,j),
     &           vol_th(1,j)
         END DO
         
         CALL DSBGV ('V','L',nth2,1,0,lap_th(:,1:nth2),2,
     &        vol_th(1,1:nth2),1,w,vec,nth2,work,info)
         
         DO j = 1, nth
            vol_th (1,j) = 
     &           (COS (theta_if (j-1)) - COS (theta_if (j))) /
     &           dtheta
         END DO
         DO j = 1, nth2 
            matrix (1:nth2,j,kk,1) = vec(1:nth2,j)
         END DO
         lambda (nth2+1:nth,kk) = DBLE (w (1:nth2))
         IF (info .NE. 0) STOP 'DSBGV failed!'
         
c     eigenfunctions (discrete associated Legendre polynomials)
         DO j = 1, nth2
            DO i = 1, nth2
               PRINT *, k, vec (i,j), 
     &              j, i, COS (theta (i)),'eigenvec (odd)'
            END DO
         END DO

c     effective l*(l+1)
         DO j = 1, nth/2
            PRINT *,'l*(l+1), odd',j, k, w (j), lambda_phi
         END DO

      END DO

      DEALLOCATE (lap_th)

      RETURN
      
      END SUBROUTINE setup_poisson_and_fft

c     ================================================================




c     ================================================================

      SUBROUTINE solve_poisson_3d (rho, phi)

c     ================================================================

      USE precision
      USE grid_cfc
      USE mo_mpi

      IMPLICIT NONE

      REAL (KIND=rk), INTENT (IN)  :: rho (nr, n_s:n_e, o_s:o_e)
      REAL (KIND=rk), INTENT (OUT) :: phi (nr, n_s:n_e, o_s:o_e)

      REAL (KIND=rk) :: x (n_s:n_e,2*nr,o_s:o_e)
      REAL (KIND=rk) :: y (n_s:n_e,2*nr,o_s:o_e,0:nth_proc/2-1)

      REAL (KIND=rk) :: src (1:nr,2)
      REAL (KIND=rk) :: scr1
      REAL (KIND=rk) :: diag (1:nr), offdiag (1:nr-1)

      DOUBLE COMPLEX :: fft1 (nr,n_s:n_e,o_s:o_e)

      INTEGER (KIND=ik) :: i, j, jj, k, kk, jk, l, l0, k0, inc, 
     &     ntmp, ndat, nth2, info
      
      REAL (KIND=rk) :: sbuf (n_s:n_e,1:2*nr,o_loc,1:nth_proc/2)
      REAL (KIND=rk) :: rbuf (n_s:n_e,1:2*nr,o_loc,1:nth_proc/2)

      DOUBLE COMPLEX :: sbufp (1:nr,n_s:n_e,o_s:o_e)
      DOUBLE COMPLEX :: rbufp (1:nr,n_s:n_e,o_s:o_e)

      INTEGER (KIND=ik) :: srce, dest, coord (2), ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      INTEGER (KIND=ik), parameter :: tag_ysndl=10000_ik
      INTEGER (KIND=ik), parameter :: tag_ysndr=10001_ik
      INTEGER (KIND=ik), parameter :: tag_zsndl=10002_ik
      INTEGER (KIND=ik), parameter :: tag_zsndr=10003_ik

      REAL (KIND=rk) :: tim1, tim2

c      CALL MPI_Barrier (cart_comm, ierr)
c      print*,'cpu_tim'
c      CALL CPU_TIME (tim1)
c      tim1 = MPI_WTime()

      nth2 = n/2

      CALL fft_forward (rho (1:nr,n_s:n_e,o_s:o_e),
     &     fft1 (1:nr,n_s:n_e,o_s:o_e), nr)

c     -----------------------------------------------------------------
c     Projection on vectors of even and odd parity

      IF (n_loc .LT. n) THEN

         l0   = (n_s-1) / n_loc
         k0   = (o_s-1) / o_loc

         coord (1) = nth_proc - 1 - l0
         coord (2) = k0
         CALL MPI_Cart_Rank (cart_comm, coord, 
     &        dest, ierr)

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
         END DO

         ndat = n_loc * nr * o_loc

         IF (l0 .LT. nth_proc/2) THEN !even parity

            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik, 
     &           cart_comm, mpistat, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              0.5d0 * (fft1 (1:nr,j,k) + rbufp (1:nr,jj,k))
            END DO

         ELSE !odd parity
            
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik,
     &           cart_comm, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              0.5d0 * (rbufp (1:nr,jj,k) - fft1 (1:nr,j,k))
            END DO
            
         END IF

      ELSE

         DO jk = 1, n_loc/2 * o_loc
            k = int((jk + n_loc/2 - 1) / (n_loc/2) )
            j = (n_s - 1) + (jk - (k - 1) * (n_loc/2))
            jj = (n_e + 1) - (jk - (k - 1) * (n_loc/2))
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
            fft1 (1:nr,j ,k) =
     &           0.5d0 * (fft1  (1:nr,jj,k) + fft1 (1:nr,j ,k))
            fft1 (1:nr,jj,k) =
     &           0.5d0 * (sbufp (1:nr,j ,k) - fft1 (1:nr,jj,k))
         END DO

      END IF

c     -----------------------------------------------------------------
c     Forward Legendre transform, local part:

      DO k = o_s, o_e

         DO i = 1, nr
            x (n_s:n_e,2*i-1,k) =
     &           DBLE  (fft1(i,n_s:n_e,k)) * vol_th(1,n_s:n_e)
            x (n_s:n_e,2*i  ,k) =
     &           DIMAG (fft1(i,n_s:n_e,k)) * vol_th(1,n_s:n_e)
         END DO

c     WARNING: this assumes that nth_proc is a least 2
         DO l = 0, nth_proc/2-1
            IF (n_e .LE. nth2) THEN
c     even parity
               CALL DGEMM ('T', 'N', n_loc, 2*nr , n_loc, 1.0d0,
     &              matrix (n_s:n_e,n_loc*l+1:n_loc*(l+1),k,0),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            ELSE
c     odd parity
               CALL DGEMM ('T', 'N', n_loc, 2*nr , n_loc, 1.0d0,
     &              matrix (n_s-nth2:n_e-nth2,
     &              n_loc*l+1:n_loc*(l+1),k,1),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            END IF
         END DO

      END DO
      
c     -----------------------------------------------------------------
c     MPI part of matrix-matrix multiplication
      
      inc = 1
      DO WHILE (inc .LE. nth_proc / 4)

         ntmp = nth_proc / (2*inc)
         l0   = (n_s-1) / (n_loc * inc)
         k0   = (o_s-1) / o_loc
         ndat = 2*nr*n_loc*o_loc*ntmp/2
            
         IF (MOD (l0, 2) .EQ. 0) THEN ! send right

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2)

            CALL MPI_Cart_shift (cart_comm, 0_ik,  inc,
     &           srce, dest, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik, 
     &           cart_comm, mpistat, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         ELSE !send left

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2)
            CALL MPI_Cart_shift (cart_comm, 0_ik,  -inc,
     &           srce, dest, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik,
     &           cart_comm, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         END IF

c     Reduction step: after the final iteration,
c     y (n_s:n_e,1:2*nr,o_s:o_e,1) contains one row
c     of the matrix-matrix-product
         DO l = 0, ntmp/2 - 1   !DO NOT PARALLELIZE THIS LOOP
            y (n_s:n_e,1:2*nr,o_s:o_e,l) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l) +
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l+1)
         END DO

         inc = 2 * inc

      END DO
      

c     -----------------------------------------------------------------
c     Now solve the tridiagonal system for each (j,k)

      DO kk = o_s, o_e

         k = index (kk-1) !frequency space index

         DO j = n_s, n_e

            DO i = 1, nr
               diag (i) = diag0 (i) - lambda (j,kk) * dv_r (i) / 
     &              r (i) ** 2
c               diag (i) = diag0 (i) - lambda (j,kk) * 
c     &              (r_if (i) - r_if (i-1))
               src (i,1) = - y (j,2*i-1,kk,0) * dv_r (i)
               src (i,2) = - y (j,2*i  ,kk,0) * dv_r (i)
            END DO

            DO i = 1, nr-1
               offdiag (i) = offdiag0 (i)
            END DO

            scr1 = SQRT (1.0d0 - 4.0d0 * lambda (j, kk))
            i = nr
            diag (i) = diag (i) + da_r (i) *
     &           (1.0d0 + scr1) / (2.0d0 * r_if (i)) *
     &           (r (i) / r_if (i)) ** scr1

            CALL DPTSV (nr, 2, diag, offdiag, src, nr, info)

            IF (info .EQ. 0) THEN
c               IF (kk .LE. nth+1-j) THEN ! filter out |m|>l
                  DO i = 1, nr
                     x (j,2*i-1,kk) = src (i,1)
                     x (j,2*i  ,kk) = src (i,2)
                  END DO
c               ELSE
c                  y (j,:) = 0.0d0
c               END IF
            ELSE
               DO i = 1, nr
                  PRINT *,i,diag(i),offdiag(i),src(i,:),lambda(j,kk)
               END DO
               STOP 'DPTSV failed!'
            END IF

         END DO !j-loop

      END DO !kk- and k-loop


c     -----------------------------------------------------------------
c     Backward Legendre transform, local part:

      DO k = o_s, o_e

         DO l = 0, nth_proc/2-1

            IF (n_e .LE. nth2) THEN
c     even parity
               CALL DGEMM ('N', 'N', n_loc, 2*nr, n_loc, 1.0d0,
     &              matrix (n_loc*l+1:n_loc*(l+1),n_s:n_e,k,0),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            ELSE
c     odd parity
               CALL DGEMM ('N', 'N', n_loc, 2*nr, n_loc, 1.0d0,
     &              matrix (n_loc*l+1:n_loc*(l+1),
     &              n_s-nth2:n_e-nth2,k,1),
     &              n_loc, x (n_s:n_e,1:2*nr,k), n_loc, 0.0d0,
     &              y (n_s:n_e,1:2*nr,k,l), n_loc)
            END IF
         END DO

      END DO

c     -----------------------------------------------------------------
c     MPI part of matrix-matrix multiplication (works exactly as before)

c     Pack
      inc = 1
      DO WHILE (inc .LE. nth_proc / 4)

         ntmp = nth_proc / (2*inc)
         l0   = (n_s-1) / (n_loc * inc)
         ndat = 2*nr*n_loc*o_loc*ntmp/2

         IF (MOD (l0, 2) .EQ. 0) THEN ! send right
c     
            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2)

            CALL MPI_Cart_shift (cart_comm, 0_ik,  inc,
     &           srce, dest, ierr)
            CALL MPI_Send (sbuf (:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik, 
     &           cart_comm, mpistat, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,1:ntmp-1:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         ELSE !send left

            sbuf (:,:,:,1:ntmp/2) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2)
            CALL MPI_Cart_shift (cart_comm, 0_ik,  -inc,
     &           srce, dest, ierr)
            CALL MPI_Recv (rbuf(:,:,:,1:ntmp/2), ndat, 
     &           MPI_DOUBLE_PRECISION, dest, 7000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbuf(:,:,:,1:ntmp/2), ndat,
     &           MPI_DOUBLE_PRECISION, dest, 7001_ik,
     &           cart_comm, ierr)
            y (n_s:n_e,1:2*nr,o_s:o_e,0:ntmp-2:2) =
     &           rbuf(:,:,:,1:ntmp/2)

         END IF

c     Reduction step: after the final iteration,
c     y (n_s:n_e,1:2*nr,o_s:o_e,1) contains one row
c     of the matrix-matrix-product
         DO l = 0, ntmp/2 - 1   !DO NOT PARALLELIZE THIS LOOP
            y (n_s:n_e,1:2*nr,o_s:o_e,l) =
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l) +
     &           y (n_s:n_e,1:2*nr,o_s:o_e,2*l+1) 
         END DO

         inc = 2 * inc

      END DO


c     -----------------------------------------------------------------
c     Pack x into complex array for inverse FFT
      DO k = o_s, o_e
         DO j = n_s, n_e
            DO i = 1, nr
               fft1 (i,j,k) = CMPLX (y (j,2*i-1,k,0), y(j,2*i,k,0),
     &              KIND(0.0d0))
            END DO
         END DO
      END DO

c     -----------------------------------------------------------------
c     Add parts with even and odd parity
      IF (n_loc .LT. n) THEN

         l0   = (n_s-1) / n_loc
         k0   = (o_s-1) / o_loc

         coord (1) = nth_proc - 1 - l0
         coord (2) = k0

         CALL MPI_Cart_Rank (cart_comm, coord, 
     &        dest, ierr)

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
         END DO

         ndat = n_loc * nr * o_loc
         
         IF (l0 .LT. nth_proc/2) THEN !even parity

            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik,
     &           cart_comm, ierr)
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik, 
     &           cart_comm, mpistat, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              fft1 (1:nr,j,k) + rbufp (1:nr,jj,k)
            END DO

         ELSE !odd parity
            
            CALL MPI_Recv (rbufp(:,:,:), ndat, 
     &           MPI_DOUBLE_COMPLEX, dest, 11000_ik, 
     &           cart_comm, mpistat, ierr)
            CALL MPI_Send (sbufp(:,:,:), ndat,
     &           MPI_DOUBLE_COMPLEX, dest, 11001_ik,
     &           cart_comm, ierr)

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               jj = (n_e + 1) - (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               fft1 (1:nr,j,k) =
     &              rbufp (1:nr,jj,k) - fft1 (1:nr,j,k)
            END DO
            
         END IF

      ELSE

         DO jk = 1, n_loc/2 * o_loc
            k = int((jk + n_loc/2 - 1) / (n_loc/2) )
            j = (n_s - 1) + (jk - (k - 1) * (n_loc/2))
            jj = (n_e + 1) - (jk - (k - 1) * (n_loc/2))
            k = k + o_s - 1
            sbufp (1:nr,j,k) = fft1 (1:nr,j,k)
            fft1 (1:nr,j ,k) =
     &           fft1  (1:nr,jj,k) + fft1 (1:nr,j ,k)
            fft1 (1:nr,jj,k) =
     &           sbufp (1:nr,j ,k) - fft1 (1:nr,jj,k)
         END DO

      END IF

c     -----------------------------------------------------------------
c     Inverse FFT
      CALL fft_inverse (fft1 (1:nr,n_s:n_e,o_s:o_e),
     &     phi(1:nr,n_s:n_e,o_s:o_e), nr)

c      tim2 = MPI_WTime()
c      CALL CPU_TIME (tim2)

c      PRINT '("Poisson: ",f8.4,x,"s")',tim2-tim1

      RETURN

      END SUBROUTINE solve_poisson_3d

c     ================================================================



#if 0

c     =================================================================

      SUBROUTINE fft_compressed (xin, yout, m, mode)

c     =================================================================

      USE precision
      USE mo_mpi

      IMPLICIT NONE

c     -----------------------------------------------------------------
c     Subroutine arguments:

      INTEGER (KIND=ik), INTENT (IN) :: m
      REAL (KIND=rk), INTENT (IN ) :: xin  (1:m,n_s:n_e,o_s:o_e)
      REAL (KIND=rk), INTENT (OUT) :: yout (1:m,n_s:n_e,o_s:o_e)

      INTEGER (KIND=ik), INTENT (in) :: mode

c     -----------------------------------------------------------------
      
      REAL (KIND=rk) :: arg

      DOUBLE COMPLEX :: y  (1:m/2,n_s:n_e,o_s-1:o_e-1)
      DOUBLE COMPLEX :: y0 (1:m/2,n_s:n_e,o_s-1:o_e-1)

      DOUBLE COMPLEX :: om, omi
c     we never have more than 2**21 cells in varphi

      INTEGER (KIND=ik) :: i, j, ji, jf, k, kk, jk, di, di0, di1,
     &     inx, nk

      INTEGER (KIND=ik) :: src, dest, ndat, ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right

c     -----------------------------------------------------------------
      
      IF (MOD (m,2) .NE. 0) STOP 'fft_compressed: m is not even!'

      ndat = (m/2) * n_loc * o_loc

c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         DO i = 1, m / 2
            y0 (i,j,k-1) = DCMPLEX (xin (2*i-1,j,k), xin(2*i,j,k))
         END DO

      END DO


      IF (mode .EQ. 1) THEN !forward transform
      
         DO kk = 1, logn
            
            di = 2 ** (logn - kk)
            nk = 2 ** (kk-1)
            arg = - DBLE (2.0d0*pi/DBLE (2*di))
            om = DCMPLX (COS (arg), SIN (arg))

            IF (2*di .LE. o_loc) THEN !no MPI communication needed

               ji = (o_s - 1) / (2*di)
               jf = (o_e - 1) / (2*di)

               DO j = ji, jf
                  di0 = 2*di*j
                  di1 = 2*di*j+di
                  DO i = 0, di-1
                     omi = DCMPLX (COS (i * arg), SIN (i * arg))
                     y (:,:,di0+i) =
     &                    y0 (:,:,di0+i) + y0 (:,:,di1+i)
                     y (:,:,di1+i) =
     &                    (y0 (:,:,di0+i) - y0 (:,:,di1+i)) * omi !om ** i
                  END DO
               END DO

            ELSE !MPI communication needed

               IF (MOD (o_s-1,di) / (di / 2) .EQ. 0) THEN

                  send_right = .TRUE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Send (y0, ndat,
     &                 MPI_DOUBLE_PRECISION, dest, 6000_ik,
     &                 cart_comm, ierr)
                  CALL MPI_Recv (y, ndat, 
     &                 MPI_DOUBLE_PRECISION, dest, 6001_ik, 
     &                 cart_comm, ierr)

                  DO i = o_s-1, o_e-1
                     inx = MOD (i, di/2)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    y0 (:,:,i) + y (:,:,i)
                  END DO

               ELSE

                  send_right = .FALSE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik, -di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Recv (y, ndat,
     &                 MPI_DOUBLE_PRECISION, dest, 6001_ik,
     &                 cart_comm, ierr)
                  CALL MPI_Send (y0, ndat,
     &                 MPI_DOUBLE_PRECISION, dest, 6000_ik, 
     &                 cart_comm, ierr)
                  DO i = o_s-1, o_e-1
                     inx = MOD (i, di/2)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    (y (:,:,i) - y0 (:,:,i)) * omi !om ** i
                  END DO
               END IF

            END IF

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 2
               y0 (:,j,k) = y (:,j,k)
            END DO

         END DO

            


c     -----------------------------------------------------------------
      ELSE IF (mode .EQ. -1) THEN ! backward transform

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 2
            y0 (:,j,k) = y (:,j,k) / DBLE (qz)
         END DO

         DO kk = logn, 1, -1

            di = 2 ** (logn - kk)
            arg = DBLE (2.0d0*pi/DBLE (2**kk))
            om = DCMPLX (COS (arg), SIN (arg))

            IF (2*di .LE. o_loc) THEN !no MPI communication needed

               ji = (o_s - 1) / (2*di)
               jf = (o_e - 1) / (2*di)

               DO j = ji, jf
                  
                  di0 = 2*di*j
                  di1 = 2*di*j+di
                  DO i = 0, di-1
                     inx = MOD(index (di0+i),2**kk)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,di0+i) = 
     &                     y0 (:,:,di0+i) + y0 (:,:,di1+i)
                     y (:,:,di1+i) =
     &                    (y0 (:,:,di0+i) - y0 (:,:,di1+i)) * omi
                  END DO
               END DO

            ELSE !MPI communication needed

               IF (MOD (o_s-1,di) / (di / 2) .EQ. 0) THEN

                  send_right = .TRUE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Send (y0, ndat,
     &                 MPI_DOUBLE_PRECISION, dest, 6000_ik,
     &                 cart_comm, ierr)
                  CALL MPI_Recv (y, ndat, 
     &                 MPI_DOUBLE_PRECISION, dest, 6001_ik, 
     &                 cart_comm, ierr)

                  DO i = o_s-1, o_e-1
                     inx = MOD(index (i),nk)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    y0 (:,:,i) + y (:,:,i)
                  END DO

               ELSE

                  send_right = .FALSE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik, -di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Recv (y, ndat,
     &                 MPI_DOUBLE_PRECISION, dest, 6001_ik,
     &                 cart_comm, ierr)
                  CALL MPI_SEND (y0, ndat, 
     &                 MPI_DOUBLE_PRECISION, dest, 6000_ik, 
     &                 cart_comm, ierr)
                  DO i = o_s-1, o_e-1
                     inx = MOD(index (i),nk)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    (y (:,:,i) - y0 (:,:,i)) * omi !om ** i
                  END DO
               END IF

            END IF

            DO jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 2
               y0 (:,j,k) = y (:,j,k)
            END DO

         END DO
            

      ELSE

         STOP 'FFT: mode not implemented!'

      END IF


c     Unpack array:
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         DO i = 1, m / 2
            yout (2*i-1,j,k) = DBLE (y (2*i-1,j,k-1))
            yout (2*i  ,j,k) = AIMAG (y (2*i  ,j,k-1))
         END DO

      END DO


      RETURN

      END SUBROUTINE fft_compressed
      
c     =================================================================

#endif      



c     =================================================================

      SUBROUTINE fft_forward (xin, yout, m)

c     =================================================================

      USE precision
      USE mo_mpi

      IMPLICIT NONE

c     -----------------------------------------------------------------
c     Subroutine arguments:

      INTEGER (KIND=ik), INTENT (IN) :: m

      REAL (KIND=rk), INTENT (IN ) :: xin  (1:m,n_s:n_e,o_s:o_e)
      DOUBLE COMPLEX, INTENT (OUT) :: yout (1:m,n_s:n_e,o_s:o_e)

c     -----------------------------------------------------------------
      
      REAL (KIND=rk) :: arg

      DOUBLE COMPLEX :: y  (1:m,n_s:n_e,o_s-1:o_e-1)
      DOUBLE COMPLEX :: y0 (1:m,n_s:n_e,o_s-1:o_e-1)

      DOUBLE COMPLEX :: om, omi

      INTEGER (KIND=ik) :: i, j, ji, jf, k, jk, di, di0, di1, inx, nk
      INTEGER (KIND=ik) :: jj, kk

      INTEGER (KIND=ik) :: src, dest, ndat, ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right

c     -----------------------------------------------------------------

      ndat = m * n_loc * o_loc

c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         y0 (1:m,j,k-1) = xin (1:m,j,k)
      END DO



      DO k = 1, logn
         
         di = 2 ** (logn - k)
         nk = 2 ** (k-1)
         arg = - DBLE (2.0d0*pi/DBLE (2*di))
         om = DCMPLX (COS (arg), SIN (arg))

         IF (2*di .LE. o_loc) THEN !no MPI communication needed
            
            ji = (o_s - 1) / (2*di)
            jf = (o_e - 1) / (2*di)
            
            DO j = ji, jf
               di0 = 2*di*j
               di1 = 2*di*j+di
               DO i = 0, di-1
                  omi = DCMPLX (COS (i * arg), SIN (i * arg))
                  y (:,:,di0+i) =
     &                 y0 (:,:,di0+i) + y0 (:,:,di1+i)
                  y (:,:,di1+i) =
     &                 (y0 (:,:,di0+i) - y0 (:,:,di1+i)) * omi
               END DO
            END DO
            
         ELSE                   !MPI communication needed

#ifdef MPI_HYDRO
            IF (MOD ((o_s-1)/di, 2) .EQ. 0) THEN
               
               send_right = .TRUE.
               CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6000_ik, cart_comm, ierr)
               CALL MPI_Recv (y,  ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6001_ik, cart_comm, mpistat, ierr) 
               
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =
     &                 y0 (:,:,i) + y (:,:,i)
               END DO
               
            ELSE
               
               send_right = .FALSE.
               CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6000_ik, cart_comm, mpistat, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6001_ik, cart_comm, ierr) 
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =
     &                 (y (:,:,i) - y0 (:,:,i)) * omi !om ** i
               END DO
            END IF
#endif /* MPI_HYDRO */

         END IF

#if 0
         DO jk = 1, n_loc * o_loc
            kk = int((jk + n_loc - 1) / n_loc )
            jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
            kk = kk + o_s - 2
            y0 (:,jj,kk-1) = y (:,jj,kk-1)
         END DO
#else
         y0 = y
#endif
         
      END DO
      
      DO jk = 1, n_loc * o_loc
         kk = int((jk + n_loc - 1) / n_loc )
         jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
         kk = kk + o_s - 1
         yout (1:m,jj,kk) = y (1:m,jj,kk-1)
      END DO


      RETURN

      END SUBROUTINE fft_forward
      
c     =================================================================



c     =================================================================

      SUBROUTINE fft_inverse (xin, yout, m)

c     =================================================================

      USE precision
      USE mo_mpi

      IMPLICIT NONE

c     -----------------------------------------------------------------
c     Subroutine arguments:

      INTEGER (KIND=ik), INTENT (IN) :: m 

      DOUBLE COMPLEX, INTENT (IN ) :: xin  (1:m,n_s:n_e,o_s:o_e)
      REAL (KIND=rk), INTENT (OUT) :: yout (1:m,n_s:n_e,o_s:o_e)

c     -----------------------------------------------------------------

      REAL (KIND=rk) :: arg

      DOUBLE COMPLEX :: y  (1:m,n_s:n_e,o_s-1:o_e-1)
      DOUBLE COMPLEX :: y0 (1:m,n_s:n_e,o_s-1:o_e-1)

      DOUBLE COMPLEX :: om, omi

      INTEGER (KIND=ik) :: i, j, ji, jf, k, jk, di, di0, di1, inx, nk
      INTEGER (KIND=ik) :: jj, kk

      INTEGER (KIND=ik) :: src, dest, ndat, ierr
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right
c     -----------------------------------------------------------------

      ndat = m * n_loc * o_loc

c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         y0 (1:m,j,k-1) = xin (1:m,j,k)
      END DO

      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2
         y0 (:,j,k) = y0 (:,j,k) / DBLE (o)
      END DO
      
      DO k = logn, 1, -1
         
         di = 2 ** (logn - k)
         nk = 2 ** (k-1)
         arg = DBLE (2.0d0*pi/DBLE (2*di))
         om = DCMPLX (COS (arg), SIN (arg))         

         IF (2*di .LE. o_loc) THEN !no MPI communication needed
#ifdef MPI_HYDRO
            ji = (o_s - 1) / (2*di)
            jf = (o_e - 1) / (2*di)
            
            DO j = ji, jf
               
               di0 = 2*di*j
               di1 = 2*di*j+di
               DO i = 0, di-1
                  omi = DCMPLX (COS (i * arg), SIN (i * arg))
                  y (:,:,di0+i) = !0.5_rk * 
     &                 (y0 (:,:,di0+i) + y0 (:,:,di1+i) * omi)
                  y (:,:,di1+i) = !0.5_rk *
     &                 (y0 (:,:,di0+i) - y0 (:,:,di1+i) * omi)
               END DO
            END DO
            
         ELSE                   !MPI communication needed

            IF (MOD ((o_s-1)/di, 2) .EQ. 0) THEN
               
               send_right = .TRUE.
               CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6000_ik, cart_comm, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              dest, 6001_ik, cart_comm, mpistat, ierr)
               
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =   !0.5_rk *
     &                 (y0 (:,:,i) + y (:,:,i) * omi)
               END DO

            ELSE

               send_right = .FALSE.
               CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
               CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6000_ik, cart_comm, mpistat, ierr)
               CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &              src, 6001_ik, cart_comm, ierr) 
               DO i = o_s-1, o_e-1
                  inx = MOD (i, di)
                  omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                  y (:,:,i) =   !0.5_rk *
     &                 (y (:,:,i) - y0 (:,:,i) * omi)
               END DO
            END IF

#endif /* MPI_HYDRO */

         END IF
         
#if 0
         DO jk = 1, n_loc * o_loc
            kk = int((jk + n_loc - 1) / n_loc )
            jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
            kk = kk + o_s - 2
            y0 (:,jj,kk) = y (:,jj,kk)
         END DO
#else
         y0 = y
#endif
      END DO


c     Unpack array:
      DO jk = 1, n_loc * o_loc
         kk = int((jk + n_loc - 1) / n_loc )
         jj = (n_s - 1) + (jk - (kk - 1) * n_loc)
         kk = kk + o_s - 1

         DO i = 1, m
            yout (i,jj,kk) = REAL (y (i,jj,kk-1), KIND=rk)
         END DO

      END DO


      RETURN

      END SUBROUTINE fft_inverse
      
c     =================================================================



c     ==================================================================

      SUBROUTINE loese_poisson (rho,phi)

c     ==================================================================
c
c     High-level driver routine for solving the Poisson equation:
c     div grad phi = 4 pi rho
c
c     Remarks:
c     - Uses multipole expansion with numerical eigenfunctions
c       of discretised Poisson equation
c     - Written for 3D case only
c     
c     rho: source density at CELL CENTRES
c     phi: solution at CELL CENTRES
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc, ONLY: equatorial_symmetry

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      USE configure

      IMPLICIT NONE

      real (kind=rk), intent (in) :: 
     &     rho (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk), intent (out) :: 
     &     phi (0:m+1,n_s-1:n_e+1,o_s-1:o_e+1)

      integer (kind=ik) :: i,j,k,jk,l

#ifdef MPI_HYDRO
      INTEGER (kind=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)

      real (kind=rk) :: sbufy (1:m, 1, o_s:o_e)
      real (kind=rk) :: rbufy (1:m, 1, o_s:o_e)
      real (kind=rk) :: sbufz (1:m, n_s-1:n_e+1, 1)
      real (kind=rk) :: rbufz (1:m, n_s-1:n_e+1, 1)

      INTEGER (kind=ik), parameter :: tag_ysndl=11000_ik
      INTEGER (kind=ik), parameter :: tag_ysndr=11001_ik
      INTEGER (kind=ik), parameter :: tag_zsndl=11002_ik
      INTEGER (kind=ik), parameter :: tag_zsndr=11003_ik
#endif


c     Call the MVM/FFT-solver
      CALL solve_poisson_3d (rho (1:m, n_s:n_e, o_s:o_e), 
     &     phi (1:m, n_s:n_e, o_s:o_e))

#ifdef MPI_HYDRO
c     Exchange first and last zone in theta and phi with the neighbouring
c     tasks

      sbufy (1:m, 1, o_s:o_e) = phi (1:m, n_s, o_s:o_e)
      CALL MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
      CALL MPI_Sendrecv (sbufy,m*o_loc,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy,
     &        m*o_loc,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndl,cart_comm,mpistat,ierr)
      phi (1:m, n_e+1, o_s:o_e) = rbufy (1:m, 1, o_s:o_e)

      sbufy (1:m, 1, o_s:o_e) = phi (1:m, n_e, o_s:o_e)
      CALL MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
      CALL MPI_Sendrecv (sbufy,m*o_loc,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndr,rbufy,
     &        m*o_loc,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndr,cart_comm,mpistat,ierr)
      phi (1:m, n_s-1, o_s:o_e) = rbufy (1:m, 1, o_s:o_e)

#endif /* MPI_HYDRO*/

c     Reflection at axis (one may be able to do better)
      IF (n_s .EQ. 1) phi (:,  0,:) = phi (:,1,:)
      IF (n_e .EQ. n) phi (:,n+1,:) = phi (:,n,:)

#ifdef MPI_HYDRO
      sbufz (1:m, n_s-1:n_e+1, 1) = phi (1:m, n_s-1:n_e+1, o_s)
      CALL MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
      IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
         CALL MPI_Send (sbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &        cart_comm,ierr)
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &        mpistat,ierr)
      ELSE
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &        mpistat,ierr)
         CALL MPI_Send (sbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &        cart_comm,ierr)
      END IF
      phi (1:m, n_s-1:n_e+1, o_e+1) = rbufz (1:m, n_s-1:n_e+1, 1)

      CALL MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
      sbufz (1:m, n_s-1:n_e+1, 1) = phi (1:m, n_s-1:n_e+1, o_e)
      IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
         CALL MPI_Send (sbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &        cart_comm,ierr)
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &        mpistat,ierr)
      ELSE
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &        mpistat,ierr)
         CALL MPI_Send (sbufz,m*(n_loc+2_ik),
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &        cart_comm,ierr)
      END IF
      phi (1:m, n_s-1:n_e+1, o_s-1) = rbufz (1:m, n_s-1:n_e+1, 1)
#else
c     Periodic boundary conditions
      phi (1:m, n_s-1:n_e+1, o_e+1) = phi (1:m, n_s-1:n_e+1, o_s)
      phi (1:m, n_s-1:n_e+1, o_s-1) = phi (1:m, n_s-1:n_e+1, o_e)
#endif /* MPI_HYDRO*/

c     Reflection at origin
      phi (0,:,:) = phi (1,:,:)
c     Extrapolate at outer boundary
      phi (m+1,:,:) = phi (m,:,:) * (r(m) / r(m+1))

      RETURN

      END SUBROUTINE loese_poisson 
      
c     ==================================================================



c     ==================================================================

      SUBROUTINE loese_vektor_poisson (s_vec, beta)

c     ==================================================================
c
c     Solves the generalized vector Poisson equation:
c     div grad beta + 1/3 grad div beta= 4 pi S
c
c     Remarks:
c     - Multipole expansion of components in global Cartesian basis
c     - Uses the Shibata scheme to treat the grad div term
c
c     s_vec: source at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c     beta:  solution at CELL CENTRES in polar spherical coordinate
c            basis (\partial_r,\partial_\theta,\partial_\varphi)
c
c     WARNING: currently restricted to axisymmetric vector field!
c
c     ------------------------------------------------------------------

      USE precision
      USE mod_poisson_aux, ONLY: beta_x, beta_y, beta_z,
     &     s_x, s_y, s_z, s_tmp, chi

#ifdef MPI_HYDRO
      USE mo_mpi
#endif
      USE size_cfc
      USE grid_cfc
      USE perm_aux_cfc, ONLY: aamp, delta_r_amp_inv,
     &     cosin_theta, sinus_theta,
     &     cosin_varphi, sinus_varphi, r_inv, rr_inv,
     &     sinus_theta_inv, delta_theta_inv,
     &     delta_varphi_inv
      USE configure

      IMPLICIT NONE

      REAL (KIND=rk), intent (in) ::
     &     s_vec (1:nr,n_s:n_e,o_s:o_e,3)
      REAL (KIND=rk), intent (out) :: 
     &     beta (0:nr+1,n_s-1:n_e+1,o_s-1:o_e+1,3)

      INTEGER (KIND=ik) :: i,j,k,jk,ji,ki,l

#ifdef MPI_HYDRO
      REAL (KIND=rk) :: sbufy (1:nr, 1, o_s:o_e, 3)
      REAL (KIND=rk) :: rbufy (1:nr, 1, o_s:o_e, 3)
      REAL (KIND=rk) :: sbufz (1:nr, n_s-1:n_e+1, 1, 3)
      REAL (KIND=rk) :: rbufz (1:nr, n_s-1:n_e+1, 1, 3)

      INTEGER (KIND=ik) :: src, dest
      INTEGER (KIND=ik) :: ierr, mpistat (MPI_STATUS_SIZE)

      integer (kind=ik), parameter :: tag_ysndl=12000_ik
      integer (kind=ik), parameter :: tag_ysndr=12001_ik
      integer (kind=ik), parameter :: tag_zsndl=12002_ik
      integer (kind=ik), parameter :: tag_zsndr=12003_ik
#endif


c Cartesian components of source vector

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      DO jk = 1, n_loc * o_loc

         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         DO i=1,nr

            s_x(i,j,k)=
     &           cosin_varphi(k)*sinus_theta(j)     *s_vec(i,j,k,1)+
     &           cosin_varphi(k)*cosin_theta(j)*r(i)*s_vec(i,j,k,2)-
     &           sinus_varphi(k)*sinus_theta(j)*r(i)*s_vec(i,j,k,3)
            s_y(i,j,k)=
     &           sinus_varphi(k)*sinus_theta(j)     *s_vec(i,j,k,1)+
     &           sinus_varphi(k)*cosin_theta(j)*r(i)*s_vec(i,j,k,2)+
     &           cosin_varphi(k)*sinus_theta(j)*r(i)*s_vec(i,j,k,3)
            s_z(i,j,k)=
     &           cosin_theta(j)     *s_vec(i,j,k,1)-
     &           sinus_theta(j)*r(i)*s_vec(i,j,k,2)
            s_tmp(i,j,k)=-r(i)*s_vec(i,j,k,1)

         END DO
      END DO

c      PRINT *,'solve_poisson_x'
      CALL loese_poisson (s_x (1:nr,n_s:n_e,o_s:o_e),
     &     beta_x (0:nr+1, n_s-1:n_e+1, o_s-1:o_e+1))
c      PRINT *,'solve_poisson_y'
      CALL loese_poisson (s_y (1:nr,n_s:n_e,o_s:o_e),
     &     beta_y (0:nr+1, n_s-1:n_e+1, o_s-1:o_e+1))
c      PRINT *,'solve_poisson_z'
      CALL loese_poisson (s_z (1:nr,n_s:n_e,o_s:o_e),
     &     beta_z (0:nr+1, n_s-1:n_e+1, o_s-1:o_e+1))
c      PRINT *,'solve_poisson_tmp'
      CALL loese_poisson (s_tmp (1:nr,n_s:n_e,o_s:o_e),
     &     chi (0:nr+1, n_s-1:n_e+1, o_s-1:o_e+1))




c     ------------------------------------------------------------------
c     Transform back to spherical polar coordinate basis
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      DO jk = 1, n_loc * o_loc

         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     treat ghost zones first
         IF (j .EQ. n_s) THEN
            ji=n_s-1
            DO i=0,nr+1
               beta(i,ji,k,1)=
     &              cosin_varphi(k)*sinus_theta(ji)*beta_x(i,ji,k)+
     &              sinus_varphi(k)*sinus_theta(ji)*beta_y(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,ji,k,2)=
     &              (cosin_varphi(k)*cosin_theta(ji)*beta_x(i,ji,k) +
     &              sinus_varphi(k)*cosin_theta(ji)*beta_y(i,ji,k) -
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,ji,k,3)=
     &              (-sinus_varphi(k) * beta_x(i,ji,k) +
     &              cosin_varphi(k) * beta_y (i,ji,k))*
     &              r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,ji,k,1)*r(i) !Shibata method
            END DO
         END IF
         IF (j .EQ. n_e) THEN
            ji=n_e+1
            DO i=0,nr+1
               beta(i,ji,k,1)=
     &              cosin_varphi(k)*sinus_theta(ji)*beta_x(i,ji,k)+
     &              sinus_varphi(k)*sinus_theta(ji)*beta_y(i,ji,k)+
     &              cosin_theta(ji)*beta_z(i,ji,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,ji,k,2)=
     &              (cosin_varphi(k)*cosin_theta(ji)*beta_x(i,ji,k) +
     &              sinus_varphi(k)*cosin_theta(ji)*beta_y(i,ji,k) -
     &              sinus_theta(ji)*beta_z(i,ji,k))*r_inv(i)
               beta(i,ji,k,3)=
     &              (-sinus_varphi(k) * beta_x(i,ji,k) +
     &              cosin_varphi(k) * beta_y (i,ji,k))*
     &              r_inv(i)*sinus_theta_inv(ji)
#endif
               chi (i,ji,k)=chi(i,ji,k)+beta(i,ji,k,1)*r(i) !Shibata method
            END DO
         END IF


         IF (k .EQ. o_s) THEN
            ki=o_s-1
            DO i=0,nr+1
               beta(i,j,ki,1)=
     &              cosin_varphi(ki)*sinus_theta(j)*beta_x(i,j,ki)+
     &              sinus_varphi(ki)*sinus_theta(j)*beta_y(i,j,ki)+
     &              cosin_theta(j)*beta_z(i,j,ki)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,j,ki,2)=
     &              (cosin_varphi(ki)*cosin_theta(j)*beta_x(i,j,ki) +
     &              sinus_varphi(ki)*cosin_theta(j)*beta_y(i,j,ki) -
     &              sinus_theta(j)*beta_z(i,j,ki))*r_inv(i)
               beta(i,j,ki,3)=
     &              (-sinus_varphi(ki) * beta_x(i,j,ki) +
     &              cosin_varphi(ki) * beta_y (i,j,ki))*
     &              r_inv(i)*sinus_theta_inv(j)
#endif
               chi (i,j,ki)=chi(i,j,ki)+beta(i,j,ki,1)*r(i) !Shibata method
            END DO
         END IF

         IF (k .EQ. o_e) THEN
            ki=o_e+1
            DO i=0,nr+1
               beta(i,j,ki,1)=
     &              cosin_varphi(ki)*sinus_theta(j)*beta_x(i,j,ki)+
     &              sinus_varphi(ki)*sinus_theta(j)*beta_y(i,j,ki)+
     &              cosin_theta(j)*beta_z(i,j,ki)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,j,ki,2)=
     &              (cosin_varphi(ki)*cosin_theta(j)*beta_x(i,j,ki) +
     &              sinus_varphi(ki)*cosin_theta(j)*beta_y(i,j,ki) -
     &              sinus_theta(j)*beta_z(i,j,ki))*r_inv(i)
               beta(i,j,ki,3)=
     &              (-sinus_varphi(ki) * beta_x(i,j,ki) +
     &              cosin_varphi(ki) * beta_y (i,j,ki))*
     &              r_inv(i)*sinus_theta_inv(j)
#endif
               chi (i,j,ki)=chi(i,j,ki)+beta(i,j,ki,1)*r(i) !Shibata method
            END DO
         END IF


c     domain interior:
         DO i=0,nr+1
               beta(i,j,k,1)=
     &              cosin_varphi(k)*sinus_theta(j)*beta_x(i,j,k)+
     &              sinus_varphi(k)*sinus_theta(j)*beta_y(i,j,k)+
     &              cosin_theta(j)*beta_z(i,j,k)
#ifndef ENFORCE_SPHERICAL_METRIC
               beta(i,j,k,2)=
     &              (cosin_varphi(k)*cosin_theta(j)*beta_x(i,j,k) +
     &              sinus_varphi(k)*cosin_theta(j)*beta_y(i,j,k) -
     &              sinus_theta(j)*beta_z(i,j,k))*r_inv(i)
               beta(i,j,k,3)=
     &              (-sinus_varphi(k) * beta_x(i,j,k) +
     &              cosin_varphi(k) * beta_y (i,j,k))*
     &              r_inv(i)*sinus_theta_inv(j)
#endif
            chi (i,j,k)=chi(i,j,k)+beta(i,j,k,1)*r(i) !Shibata method
         END DO

c     ------------------------------------------------------------------
c     Fill ghost zones (auxiliary potential chi only):
         chi(0,j,k)=chi(1,j,k)
         IF (j .EQ. 1) THEN
            chi(:,j-1,k)=chi(:,j,k)
         END IF
         IF (j .EQ. n) THEN
            chi(:,j+1,k)=chi(:,j,k)
         END IF

      END DO


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      DO jk = 1, n_loc * o_loc

         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     ------------------------------------------------------------------
c     Shibata method: beta = W + grad (chi + r * W)

         DO i=1,nr
            beta(i,j,k,1)=beta(i,j,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(i+1,j,k)-aamp(i)*chi(i-1,j,k)+
     &           (aamp(i)-1.0_rk)*chi(i,j,k))
         END DO
         
         beta(0,j,k,1)=beta(0,j,k,1)-0.125_rk*
     &           delta_r_amp_inv(i)*
     &           (chi(2,j,k)-aamp(1)*chi(1-1,j,k)+
     &           (aamp(1)-1.0_rk)*chi(1,j,k))
         beta(nr+1,j,k,1)=beta(nr+1,j,k,1)-0.125_rk*
     &        (chi(nr+1,j,k)-chi(nr,j,k))/delta_r_if(nr)

#ifndef ENFORCE_SPHERICAL_METRIC
         DO i=1,nr+1
            beta(i,j,k,2)=beta(i,j,k,2)-0.0625_rk*
     &           delta_theta_inv*rr_inv(i)*
     &           (chi(i,j+1,k)-chi(i,j-1,k))
         END DO
         beta(0,j,k,2)=beta(0,j,k,2)+0.0625_rk* !dchi/dtheta=0 on the rotation axis
     &        delta_theta_inv*rr_inv(i)*
     &        (chi(1,j+1,k)-chi(1,j-1,k))

         DO i=1,nr+1
            beta(i,j,k,3)=beta(i,j,k,2)-0.0625_rk*
     &           delta_varphi_inv*sinus_theta_inv(j)*rr_inv(i)*
     &           (chi(i,j,k+1)-chi(i,j,k-1))
         END DO
         beta(0,j,k,2)=beta(0,j,k,2)+0.0625_rk* !dchi/dtheta=0 on the rotation axis
     &        delta_varphi_inv*sinus_theta_inv(j)*rr_inv(i)*
     &        (chi(1,j,k+1)-chi(1,j,k-1))

#endif
      END DO


#ifdef MPI_HYDRO
c     Exchange first and last zone in theta and phi with the neighbouring
c     tasks

      sbufy (1:nr, 1, o_s:o_e,1:3) = beta (1:nr, n_s, o_s:o_e,1:3)
      CALL MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
      CALL MPI_Sendrecv (sbufy,nr*o_loc*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy,
     &        nr*o_loc*3_ik,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndl,cart_comm,mpistat,ierr)
      beta (1:nr, n_e+1, o_s:o_e, 1:3) = rbufy (1:nr, 1, o_s:o_e, 1:3)

      sbufy (1:nr, 1, o_s:o_e, 1:3) = beta (1:nr, n_e, o_s:o_e, 1:3)
      CALL MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
      CALL MPI_Sendrecv (sbufy,nr*o_loc*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_ysndr,rbufy,
     &        nr*o_loc*3_ik,MPI_DOUBLE_PRECISION,
     &        src,tag_ysndr,cart_comm,mpistat,ierr)
      beta (1:nr, n_s-1, o_s:o_e, 1:3) = rbufy (1:nr, 1, o_s:o_e, 1:3)

#endif /* MPI_HYDRO*/

c     Reflection at axis (one may be able to do better)
      IF (n_s .EQ. 1) THEN
         beta (:,  0,:,1) = beta (:,1,:,1)
         beta (:,  0,:,2) =-beta (:,1,:,2)
         beta (:,  0,:,3) = beta (:,1,:,3)
      END IF
      IF (n_e .EQ. n) THEN
         beta (:,n+1,:,1) = beta (:,n,:,1)
         beta (:,n+1,:,2) =-beta (:,n,:,2)
         beta (:,n+1,:,3) = beta (:,n,:,3)
      END IF

#ifdef MPI_HYDRO
      sbufz (1:nr,n_s-1:n_e+1,1,1:3) = beta (1:nr,n_s-1:n_e+1, o_s,1:3)
      CALL MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
      IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
         CALL MPI_Send (sbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &        cart_comm,ierr)
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &        mpistat,ierr)
      ELSE
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,src,tag_zsndl,cart_comm,
     &        mpistat,ierr)
         CALL MPI_Send (sbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndl,
     &        cart_comm,ierr)
      END IF
      beta (1:nr,n_s-1:n_e+1,o_e+1,1:3) = rbufz (1:nr,n_s-1:n_e+1,1,1:3)

      sbufz (1:nr, n_s-1:n_e+1,1,1:3) = beta (1:nr,n_s-1:n_e+1,o_e,1:3)
      CALL MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
      IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
         CALL MPI_Send (sbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &        cart_comm,ierr)
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &        mpistat,ierr)
      ELSE
         CALL MPI_Recv (rbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &        mpistat,ierr)
         CALL MPI_Send (sbufz,m*(n_loc+2_ik)*3_ik,
     &        MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &        cart_comm,ierr)
      END IF
      beta (1:nr, n_s-1:n_e+1,o_s-1,1:3)= rbufz (1:nr,n_s-1:n_e+1,1,1:3)
#else
c     Periodic boundary conditions
      beta (1:nr, n_s-1:n_e+1, o_e+1) = beta (1:nr, n_s-1:n_e+1, o_s)
      beta (1:nr, n_s-1:n_e+1, o_s-1) = beta (1:nr, n_s-1:n_e+1, o_e)
#endif /* MPI_HYDRO*/

c     Reflection at origin
      beta (0,:,:,1) = - beta (1,:,:,1)
      beta (0,:,:,2) = - beta (1,:,:,2)
      beta (0,:,:,3) =   beta (1,:,:,3)
c     Extrapolate at outer boundary
      beta (nr+1,:,:,:) = beta (nr,:,:,:) * (r(nr) / r(nr+1))


      RETURN

      END SUBROUTINE loese_vektor_poisson

c     ==================================================================



      END MODULE poisson3d

c     ================================================================
#endif /* (DIMN==3) && !defined(ENFORCE_SPHERICAL_METRIC) */
