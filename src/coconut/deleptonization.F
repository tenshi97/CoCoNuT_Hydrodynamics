c     ==================================================================

      MODULE delep_scheme

c     ==================================================================

      USE precision

      IMPLICIT NONE

      SAVE

      PRIVATE
      PUBLIC :: read_in_y_e_profile, deleptonization


c     ------------------------------------------------------------------
c     Parameters:
      real (kind=rk), parameter :: rho_trapping = 2.0e12_rk
      real (kind=rk), parameter :: e_nu_esc = 10.0_rk

c     ------------------------------------------------------------------
c     Tabulated data:

      real (kind=rk), allocatable :: rho_data(:), y_e_data(:)
      integer (kind=ik) :: m_y_e_data


      CONTAINS


c     ==================================================================

      SUBROUTINE read_in_y_e_profile

c     ==================================================================
      
      USE precision

      USE size_cfc
      USE parameters_cfc
      
      IMPLICIT NONE

      integer (kind=ik) :: i, istat
      
      print *, 'Reading Y_e(rho) profile.'
      
 10   format (1x, e16.6, e16.6)
      
      open (unit = 10, file = 'tables/y_e_profile.dat',
     &     status = 'unknown', access = 'sequential',
     &     form = 'formatted')

      m_y_e_data = 366!595

      allocate (rho_data (1:m_y_e_data), y_e_data (1:m_y_e_data),
     &     stat = istat)

      do i = 1, m_y_e_data
         read (10, 10) rho_data (i), y_e_data (i)
         print *,i,rho_data(i),y_e_data(i)
      enddo

      print *, 'Data points read:',m_y_e_data
      
      close (10)

      return

      END SUBROUTINE read_in_y_e_profile

c     ==================================================================




c     ==================================================================

      SUBROUTINE deleptonization (bounce)

c     ==================================================================
      
      USE precision

      USE size_cfc

      USE conserved_cfc
      USE hydro_primitives_cfc
      USE metric_cfc, ONLY: alpha, g_up_11, g_up_22, g_up_33, sqrt_gamma
      USE parameters_cfc
      USE perm_aux_cfc, ONLY: aamp, delta_r_amp_inv
      USE phycon
      USE mo_mpi
      USE eos_sn2
      USE cons_check
      USE gr_initial_setup

      USE configure

      IMPLICIT NONE
      real(kind=rk)     :: eos_self(2), eos_children(2)
      integer (kind=ik) ::  i, j, k, jk, i_nuc
      
      real (kind=rk) :: mu_nu_with_offset, eta, p_nu_factor
      real (kind=rk) :: delta_y_e, delta_entropy
      real (kind=rk) :: p_nu (0 : m + 1)
      real (kind=rk) :: rho_h_ww
      real (kind=rk) :: delta_v_r
      real (kind=rk) :: alpha_trap, delta_v_r_trap

      real (kind=rk) :: den_wrk(m),ei_wrk(m),tem_wrk(m),pre_wrk(m),
     &     s_wrk(m),gam_wrk(m),xnu_wrk(m,config%qn)

! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=rk) ::  ccu(m), cce(m), ccn(m), ccp(m)
      real (kind=rk) ::  dmy1(m), dmy2(m,2) !dummys
! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical, intent (in) :: bounce
      logical :: bounce2, ler
#ifdef MPI_HYDRO
      integer :: ierr
#endif

      p_nu_factor = 4.0_rk * pi * pc_kb ** 4 /
     &     (3.0_rk * (2.0_rk * pi * wc_hc * pc_meverg) ** 3) *
     &     pc_egeo

#if 0
      bounce = any(rho(1:m,n_s:n_e,o_s:o_e) .gt. rho_trapping*pc_ggeo
     &     .and. entropy(1:m,n_s:n_e,o_s:o_e) .gt. 3.0_rk)
      bounce = bounce .and. any(rho(1:m,n_s:n_e,o_s:o_e) .gt. 2.7e14_rk*
     &     pc_ggeo)

c      bounce = any(rho(1:m,n_s:n_e,o_s:o_e) .gt. rho_trapping*pc_ggeo
c     &     .and. entropy(1:m,n_s:n_e,o_s:o_e) .gt. 3.0_rk)
c      bounce = any(rho(1:m,n_s:n_e,o_s:o_e) .gt. 1.0e14_rk)

#ifdef MPI_HYDRO
      call MPI_Allreduce (bounce, bounce2, 1_ik,
     &     MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      bounce = bounce2
#endif
#endif
                     
#if defined(OPEN_MP_2D) || defined(OPEN_MP_3D)
!$OMP PARALLEL DO
!$OMP& SHARED (p_nu_factor,bounce,n_s,n_e,o_s,o_e,n_loc,o_loc,pi,
!$OMP& rho,v_1,v_2,v_3,v_squared,w,eps,h,t,entropy,cpot,delta_t,
!$OMP& p,gamm,c_sound_squared,
!$OMP& xnnu,alpha,g_up_11,g_up_22,g_up_33,aamp,delta_r_amp_inv,config,
!$OMP& d_cap_hat,tau_hat,s_1_hat,s_2_hat,s_3_hat,d_cap_xnu_hat,
#ifdef CFC_MHD
!$OMP& b_1,b_2,b_3,psi,
#endif      
!$OMP& m)
!$OMP& PRIVATE(i,j,k,jk,mu_nu_with_offset,eta,delta_y_e,
!$OMP& delta_entropy,p_nu,delta_v_r,den_wrk,ei_wrk,
!$OMP& tem_wrk,pre_wrk,s_wrk,gam_wrk,xnu_wrk,ccu,cce,ccp,ccn,
!$OMP& dmy1,dmy2,alpha_trap,delta_v_r_trap,ler, eos_self,eos_children)
!$OMP& DEFAULT(NONE)
#endif
      do jk = 1, n_loc * o_loc         
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c         if (.not.bounce) then

            den_wrk(1:m) = rho(1:m,j,k)*pc_geog 
            ei_wrk (1:m) = (eps(1:m,j,k)*rho(1:m,j,k))*pc_geoe
            tem_wrk(1:m) = t      (1:m,j,k)
            s_wrk  (1:m) = entropy(1:m,j,k)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
            xnu_wrk(1:m,1:config%qn)=xnnu (1:m,j,k,1:config%qn)
            
            do i = 1, m
               
               call calculate_delta_y_e
     &              (den_wrk (i), xnu_wrk (i,config%qn), delta_y_e)

               delta_y_e = min (delta_y_e, 0.0_rk)
               xnu_wrk (i, config%qn) = xnu_wrk (i, config%qn) +
     &              delta_y_e
               
               mu_nu_with_offset = cpot (i, j, k, 1) - e_nu_esc
               mu_nu_with_offset = max(mu_nu_with_offset,0.0_rk)
               
               if (den_wrk (i).lt.rho_trapping.and.
     &              mu_nu_with_offset.gt.0.0_rk) then
#if 1                        
                  delta_entropy =
     &                 - delta_y_e * mu_nu_with_offset /
     &                 (t (i, j, k) * pc_kmev)
                  
                  s_wrk (i) = s_wrk (i) + delta_entropy 
#else
                  ei_wrk(i) = ei_wrk(i) +
     &                 min(e_nu_esc-cpot(i,j,k,1),0.0_rk) * pc_meverg *
     &                 delta_y_e * den_wrk(i) / pc_mb
#endif
               end if
            enddo

#if 1
            call eos(den_wrk(1:m),tem_wrk(1:m),xnu_wrk(1:m,:),
     &           dmy1(1:m),dmy2(1:m,:),ei_wrk(1:m),pre_wrk(1:m),
     &           gam_wrk(1:m),s_wrk(1:m),ccu(1:m),cce(1:m),
     &           ccn(1:m),ccp(1:m),eos_self, eos_children, 
     &           mode=4,nsemode=0,ler=ler)
#else
            call eos(den_wrk(1:m),tem_wrk(1:m),xnu_wrk(1:m,:),
     &           dmy1(1:m),dmy2(1:m,:),ei_wrk(1:m),pre_wrk(1:m),
     &           gam_wrk(1:m),s_wrk(1:m),ccu(1:m),cce(1:m),
     &           ccn(1:m),ccp(1:m),eos_self, eos_children,
     &           mode=3,nsemode=0,ler=ler)
#endif

            p       (1:m,j,k)  = pre_wrk(1:m)*pc_egeo
            t       (1:m,j,k)  = tem_wrk(1:m)
            eps     (1:m,j,k)  = ei_wrk (1:m)/rho(1:m,j,k)*pc_egeo
            entropy (1:m,j,k)  = s_wrk (1:m)
            gamm    (1:m,j,k)  = gam_wrk(1:m)
            xnnu(1:m,j,k,1:config%qn) = xnu_wrk(1:m,1:config%qn)
            cpot (1:m,j,k,1)   = ccu(1:m)
            cpot (1:m,j,k,2)   = cce(1:m)
            cpot (1:m,j,k,3)   = ccn(1:m)
            cpot (1:m,j,k,4)   = ccp(1:m)
#ifdef CFC_TRANSPORT2
            h    (1:m,j,k)     = 1.0_rk+
     &           eps(1:m,j,k)+p(1:m,j,k)/rho(1:m,j,k)
            c_sound_squared(1:m,j,k) =
     &           gam_wrk(1:m)*p(1:m,j,k)/ (rho(1:m,j,k)*h(1:m,j,k))
#else
            c_sound_squared(1:m,j,k) =
     &           gam_wrk(1:m)*p(1:m,j,k)/rho(1:m,j,k)
#endif /* CFC_TRANSPORT2 */
            
            if (ler) stop 'ERROR in eos (deleptonization)'
c         end if
               
         do i = 1, m
               
            eta = cpot (i, j, k, 1) / (t (i, j, k) * pc_kmev) 
            
c            p_nu (i) = p_nu_factor *
c     &              t (i, j, k) ** 4 *
c     &           (0.5_rk * eta ** 2 *
c     &           (pi ** 2 + 0.5_rk * eta ** 2) +
c     &           7.0_rk / 60.0_rk * pi ** 4)

            p_nu (i) = p_nu_factor *
     &              t (i, j, k) ** 4 *
     &           (0.5_rk * eta ** 2 *
     &           (pi ** 2 + 0.5_rk * eta ** 2) +
     &           7.0_rk / 60.0_rk * pi ** 4 *
     &           (1.0_rk - 0.5_rk * exp (-abs(eta))))


         end do

         p_nu (0) = p_nu (1)
         p_nu (m + 1) = p_nu (m)
         

         delta_v_r_trap = 0.0_rk
         alpha_trap = 1.0_rk
#ifdef CFC_TRANSPORT2
c     GR version
#if 1
         do i = 1, m
               
            if (den_wrk (i) .gt. rho_trapping) then
               delta_v_r = 
     &              - delta_r_amp_inv (i) *
     &              (p_nu (i + 1) +
     &                 (aamp (i) - 1.0_rk) * p_nu (i) -
     &              aamp (i) * p_nu (i - 1)) /
     &              (rho (i, j, k) * h (i, j, k)) *
     &              alpha (i, j, k)
               delta_v_r_trap = delta_v_r
               alpha_trap = alpha (i, j, k)
            else
               delta_v_r = delta_v_r_trap * 
     &              (alpha (i, j, k) / alpha_trap) *
     &              (den_wrk (i) / rho_trapping)
            end if

            v_1 (i, j, k) = v_1 (i, j, k) +
     &           delta_t * delta_v_r  
            
            v_squared (i, j, k) =
     &           g_up_11 (i, j, k) * v_1 (i, j, k) ** 2 +
     &           g_up_22 (i, j, k) * v_2 (i, j, k) ** 2 +
     &           g_up_33 (i, j, k) * v_3 (i, j, k) ** 2

            w (i, j, k) =
     &           1.0_rk / sqrt (1.0_rk - v_squared (i, j, k))

c            s_1_hat (i, j, k) = s_1_hat (i, j, k) + delta_t *
c     &           rho (i, j, k) * phi (i, j, k) ** 2 *
c     &           delta_v_r_spherical_dt (i) *
c     &           sqrt_gamma (i, j, k)
               
c            tau_hat (i, j, k) = tau_hat (i, j, k) + delta_t *
c     &           rho (i, j, k) * phi (i, j, k) ** 2 *
c     &           v_1 (i, j, k) * g_up_11 (i, j, k) *
c     &           delta_v_r_spherical_dt (i) *
c     &           sqrt_gamma (i, j, k)

            rho_h_ww = rho (i, j, k) * h (i, j, k) *
     &           w (i, j, k) ** 2
            
            d_cap_hat (i, j, k) =  rho (i, j, k) * w (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_1_hat   (i, j, k) =  rho_h_ww * v_1 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_2_hat   (i, j, k) =  rho_h_ww * v_2 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_3_hat   (i, j, k) =  rho_h_ww * v_3 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            tau_hat   (i, j, k) = (rho_h_ww - p (i, j, k)) *
     &           sqrt_gamma (i, j, k) - d_cap_hat (i, j, k)

         end do

         do i_nuc = 1, config%qn
            do i = 1, m
               d_cap_xnu_hat (i, j, k, i_nuc) = d_cap_hat (i, j, k) *
     &              xnnu (i, j, k, i_nuc)
            end do
         enddo
#endif
#else /* CFC_TRANSPORT2 */
c     Newtonian version
#if 1
         do i = 1, m
                       
            if (den_wrk (i) .gt. rho_trapping) then
               delta_v_r = 
     &              - delta_r_amp_inv (i) *
     &              (p_nu (i + 1) +
     &                 (aamp (i) - 1.0_rk) * p_nu (i) -
     &              aamp (i) * p_nu (i - 1)) /
     &              rho (i, j, k)
               delta_v_r_trap = delta_v_r
               alpha_trap = alpha (i, j, k)
            else
               delta_v_r = delta_v_r_trap * 
     &              (den_wrk (i) / rho_trapping)
            end if

            v_1 (i, j, k) = v_1 (i, j, k) +
     &           delta_t * delta_v_r  
            
            v_squared (i, j, k) =
     &           g_up_11 (i, j, k) * v_1 (i, j, k) ** 2 +
     &           g_up_22 (i, j, k) * v_2 (i, j, k) ** 2 +
     &           g_up_33 (i, j, k) * v_3 (i, j, k) ** 2

            w (i, j, k) = 1.0_rk

            h    (i,j,k)     = eps(i,j,k) + p(i,j,k) / rho(i,j,k) +
     &        0.5_rk * v_squared (i,j,k)

            d_cap_hat (i, j, k) = rho (i, j, k)
            s_1_hat   (i, j, k) = rho (i, j, k) * v_1 (i, j, k)
            s_2_hat   (i, j, k) = rho (i, j, k) * v_2 (i, j, k)
            s_3_hat   (i, j, k) = rho (i, j, k) * v_3 (i, j, k)
            tau_hat   (i, j, k) = rho (i, j, k) * h (i, j, k) -
     &           p (i, j, k)
#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
     &           + 0.5_rk * (b_1 (i, j, k) ** 2 + b_2 (i, j, k) ** 2 +
     &           b_3 (i, j, k) ** 2)
#ifndef CONSTRAINED
     &           + 0.5_rk * psi (i, j, k) ** 2
#endif
#endif /* CFC_MHD and not PASSIVE_MHD */

         end do

         do i_nuc = 1, config%qn
            do i = 1, m
               d_cap_xnu_hat (i, j, k, i_nuc) = d_cap_hat (i, j, k) *
     &              xnnu (i, j, k, i_nuc)
            end do
         enddo
#endif
#endif /* CFC_TRANSPORT2 */
      end do                    !j,k-loop

c      call initialize_conserved_quantities
      
      return

      END SUBROUTINE deleptonization

c     ==================================================================





c     ==================================================================

      SUBROUTINE calculate_delta_y_e (rho_tmp, y_e_tmp, delta_y_e)

c     ==================================================================
      
      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc

      IMPLICIT NONE
      
      integer (kind=ik) :: ii
      
      real (kind=rk), intent (in) ::  rho_tmp, y_e_tmp
      real (kind=rk), intent (out) :: delta_y_e
      
      real (kind=rk)::  log_rho, x_rho, y_e_from_profile

      delta_y_e=0.0_rk

      do ii = 1, m_y_e_data - 1
         if (rho_data (ii).le.rho_tmp
     &        .and.rho_tmp.lt.rho_data (ii + 1)) then
            y_e_from_profile =
     &           y_e_data (ii) + 
     &           (rho_tmp - rho_data (ii)) /
     &           (rho_data (ii + 1) - rho_data (ii)) *
     &           (y_e_data (ii + 1) - y_e_data (ii))
            delta_y_e = min (0.0_rk,
     &           y_e_from_profile - y_e_tmp)
         endif
      enddo
            
      

      return

      END SUBROUTINE calculate_delta_y_e

c     ==================================================================





      END MODULE delep_scheme

c     ==================================================================
