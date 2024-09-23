#define ANALYTIC_SCATTERING

 module rice_phase
  implicit none

contains
  subroutine phase_step(f, dfdt_adv, f_eq, kappa_a, kappa_s, dt, diffusive_limit)
    use rice_config,       only: npsi, nmu, neps, nflav, clight
    use rice_grid,         only: domega_scaled
#ifndef ANALYTIC_SCATTERING
    use rice_matrixsolve,  only: solve_matrix
#endif

    real, intent(inout) :: f       (1:nflav,1:neps,1:npsi,-nmu:nmu)
    real, intent(in)    :: dfdt_adv(1:nflav,1:neps,1:npsi,-nmu:nmu)

    real, intent(in)    :: f_eq   (1:nflav,1:neps)
    real, intent(in)    :: kappa_a(1:nflav,1:neps)
    real, intent(in)    :: kappa_s(1:nflav,1:neps)
    real, intent(in)    :: dt
    logical, optional, intent(in) :: diffusive_limit

    real :: rhs(1:nflav,1:neps,npsi*(1+2*nmu))

#ifdef ANALYTIC_SCATTERING
    real :: ka_dt(1:nflav,1:neps)
    real :: ks_dt(1:nflav,1:neps)
    real :: kt_dt(1:nflav,1:neps)
    real :: scr1(1:nflav,1:neps)
    real :: denomi(1:nflav,1:neps)
#else
    integer :: index2, m, n
    real    :: matrix(npsi*(1+2*nmu),npsi*(1+2*nmu))
#endif

    integer :: j, k
    integer :: index

    logical :: diffusion

    scr1 = 0.0

    diffusion = .false.
    if (present(diffusive_limit)) diffusion = diffusive_limit

    if (diffusion) then ! Assume diffusion limit

      do j = -nmu, nmu
        do k = 1, npsi
           index = k + (j+nmu)*npsi
           rhs(:,:,index) = kappa_a(:,:) * clight * f_eq(:,:) + dfdt_adv(:,:,k,j)
           scr1(:,:) = scr1(:,:) + rhs(:,:,index) * domega_scaled(k,j)
        enddo
      enddo

      denomi(:,:) = 1.0 / (clight*(kappa_a(:,:) + kappa_s(:,:)))

      do j = -nmu, nmu
       do k = 1, npsi
          index = k + (j+nmu)*npsi
          f(:,:,k,j) = ( rhs(:,:,index) + kappa_s(:,:)*clight*denomi(:,:) ) * denomi(:,:)
       enddo
      enddo

    else ! Normal implicit solve

      do j = -nmu, nmu
        do k = 1, npsi
          index = k + (j+nmu)*npsi
          rhs(:,:,index) = f(:,:,k,j) + dt * (kappa_a(:,:) * clight * f_eq(:,:) + dfdt_adv(:,:,k,j))
          scr1(:,:) = scr1(:,:) + rhs(:,:,index) * domega_scaled(k,j)
        enddo
      enddo

      ka_dt(:,:) = kappa_a(:,:) * clight * dt
      ks_dt(:,:) = kappa_s(:,:) * clight * dt
      kt_dt(:,:) = ka_dt(:,:) + ks_dt(:,:)
      denomi(:,:) = 1.0 / ((1.0 + ka_dt(:,:)) * (1.0 + kt_dt(:,:)))

      do j = -nmu, nmu
       do k = 1, npsi
          index = k + (j+nmu)*npsi
          f(:,:,k,j) = rhs(:,:,index) / (1.0 + kt_dt(:,:)) + ks_dt(:,:) * scr1(:,:) * denomi(:,:)
       enddo
      enddo

    endif


! #else /* ANALYTIC_SCATTERING */
!         do j = -nmu, nmu
!            do k = 1, npsi
!
!               index = k + (j+nmu)*npsi
!               rhs(index) = f(l,k,j,i) + dt * (kappa_a(l,i) * clight * f_eq(l,i) + dfdt_adv(l,k,j,i))
!
!             do m = -nmu, nmu
!               do n = 1, npsi
!                 index2 = n + (m+nmu)*npsi
!                 if (index2 == index) then
!                   matrix(index, index2) = 1.0 + dt * clight * (kappa_a(l,i) + kappa_s(l,i) * (1.0 - domega_scaled(n,m)))
!                 else
!                   matrix(index, index2) = dt * (- kappa_s(l,i) * clight * domega_scaled(n,m))
!                 endif
!               enddo
!             enddo
!
!           enddo
!         enddo
!
!         call solve_matrix(matrix_size, matrix, rhs)
!         ! rhs becomes lhs
!
!         do j = -nmu, nmu
!           index = (j+nmu)*npsi
!           f(l,:,j,i) = rhs(index+1:index+npsi)
!         enddo
! #endif /* ANALYTIC_SCATTERING */

  end subroutine phase_step

end module rice_phase
