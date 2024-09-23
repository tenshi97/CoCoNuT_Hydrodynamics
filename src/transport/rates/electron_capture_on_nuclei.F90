#undef IBM_OPT

module lms_rates

implicit none

private

public abso_nuci_lms

contains

function abso_nuci_lms(em, isi, isf, hydro_quantities) result(aab)


! Purpose: calculate opacities for
!          absorption of neutrinos on heavy nuclei
!          according to an NSE table provided by K.Langanke, J.M.Sampaio
!
! Reference: Langanke, Martinez-Pinedo, Sampaio, PRC~64, 055801 (2001)
!
! This subroutine provides the rates that must be multiplied with
! the mass-fraction of the heavy nucleus. Per default this is
! taken from the LMS-NSE-table but it can be switched to the appropiate
! EoS
!
! For convienience as much as possible is identcal to the original
! subroutine abso_nuci_lms
!
! Input:
!        em  :        Neutrino energies
!        aab :        absorption opacity (due to other interactions)
!        den :        baryon density [g/cm^3]
!        tem :        temperature [MeV]
!        yee :        electron fraction
!        cel :        chem potential of electrons (including restmass)
!        cpr :        chem potential of protons   (without restmass)
!        cne :        chem potential of neutrons  (without restmass)
!        q_ga:        quenching factor \in [0,1] of g_A
!
!  only required for zones outside the LMS table, where we switch to
!   a (modified) FFN description:
!        dnuc:        NUMBER density of nuclei AZ,AN,AA
!        az,aa:       Z,A
!
! Output:
!        aab:         aab + opacities for abs. of neutrinos on nuclei
!
! Author: A.Marek (MPA) & M.Rampp (MPA)

  use precision

  use phycon
  use machcons
!  use neutrinotypes
  use rates_data_type, only : lms_rates_table
  use data_lms_rates, only : load_lmstable, iye_tb
  use interpol_eul, only : sort_vec
  use abort
  use timings
  use configure
  use mo_mpi
  use configure
  use mod_rate_types, only : hydro_quantities_t, absorption_opacity_t
  implicit none


! LOCAL variables that are not in modules
  integer(kind=ik), intent(in)          :: isi, isf

  type(hydro_quantities_t), intent(in)   :: hydro_quantities(:)
  real(kind=rk), intent(in)              :: em(:)
  real(kind=rk), dimension(size(hydro_quantities), &
                size(em))                :: aab
!      type(absorption_opacity_t),       &
!      dimension(size(hydro_quantities), &
!                size(em))                :: opacity

  real(kind=rk), dimension(size(hydro_quantities)) :: beta, qeff, chel, chep, chen, &
                                         yhvy, rtot, rffn, dist

  real(kind=rk), dimension(size(hydro_quantities)) :: log_den, chenu
  real(kind=rk), &
  dimension(size(hydro_quantities),size(em))       :: spectra

  integer(kind=ik), &
  dimension(size(hydro_quantities))                :: indx_dn, indx_tm, indx_yel, indx_yer

  real(kind=rk), dimension(size(hydro_quantities)) :: in_lmstb, in_lmstb_den,   &
                                         in_lmstb_tem, in_lmstb_ye


  real(kind=rk)                       :: scrtch(iye_tb)

  logical, save                       :: first=.true.
  integer(kind=ik)                    :: idn, idn1, itm, itm1, iyel,       &
                                         iyel1, iyer, iyer1, nnmax, nemax, &
                                         i, ie
  real(kind=rk) ::                       earge, eargn, th, x_dn, x_tm,     &
                                         x_yel, x_yer, cnu_m,              &
                                         dstd, dstt, dsty, rffn_e, rlms_e, &
                                         qffn, fnp, fnh

  real(kind=rk)                       :: x_ld

  real(kind=rk), parameter            :: delta=3._rk    !FFN Excitation energy [MeV]
  real(kind=rk), parameter            :: den_bnd=1.e8_rk, tem_bnd=0.05_rk, yee_bnd=0.01_rk
  real(kind=rk), parameter            :: cnst_phsp=2._rk*pc_pi**2*wc_hc**3/pc_cl



#ifdef HPM
#error "where is this include?"
!#include "./Include/f_hpm.h"
#endif

!  call timer%start("abso_nuci_lms")


  if (first) then

     call load_lmstable

     first=.false.

     if (config%use_low_den_electron_captures) then
        print *," "
        print *,"Attention the electron captures are switched on "
        print *,"in the low-density regime. Hopefully you know "
        print *,"what you are doing!!!"
        print *," "
     endif
  endif

  aab(:,:) = 0.0_rk

  nnmax=size(hydro_quantities)
  nemax=size(em)

  beta(:)=1._rk/hydro_quantities(:)%tem

  if (config%lms_spectra) then
     log_den(:) = log(hydro_quantities(:)%rho)
  endif

  where(hydro_quantities(:)%rho .ge. lms_rates_table(1)%romin .and. &
        hydro_quantities(:)%rho .le. lms_rates_table(1)%romax .and. &
        hydro_quantities(:)%tem .ge. lms_rates_table(1)%ttmin .and. &
        hydro_quantities(:)%tem .le. lms_rates_table(1)%ttmax .and. &
        hydro_quantities(:)%ye  .ge. lms_rates_table(1)%yemin .and. &
        hydro_quantities(:)%ye  .le. lms_rates_table(1)%yemax)
     in_lmstb=  1._rk
  elsewhere
     in_lmstb=  -1._rk
  endwhere

  where(hydro_quantities(:)%rho .ge. lms_rates_table(1)%romin .and. &
        hydro_quantities(:)%rho .le. lms_rates_table(1)%romax)
     in_lmstb_den=1._rk
  elsewhere
     in_lmstb_den=-1._rk
  endwhere

  where(hydro_quantities(:)%tem .ge. lms_rates_table(1)%ttmin .and. &
        hydro_quantities(:)%tem .le. lms_rates_table(1)%ttmax)
     in_lmstb_tem =1._rk
  elsewhere
     in_lmstb_tem=-1._rk
  endwhere

  where(hydro_quantities(:)%ye .ge. lms_rates_table(1)%yemin .and. &
        hydro_quantities(:)%ye .le. lms_rates_table(1)%yemax)
     in_lmstb_ye = 1._rk
  elsewhere
     in_lmstb_ye = -1._rk
  endwhere

! search for indices in the table (note: *_tb is defined real(kind=rk))

  call sort_vec(lms_rates_table(1)%rho,lms_rates_table(1)%nro,1,hydro_quantities%rho,indx_dn,nnmax,1)
  call sort_vec(lms_rates_table(1)%tem,lms_rates_table(1)%ntt,1,hydro_quantities%tem,indx_tm,nnmax,1)

  do i=1,nnmax
     idn =max(indx_dn(i) ,1) ; idn1 =min(idn+1,lms_rates_table(1)%nro)
     scrtch=lms_rates_table(1)%ye(:,idn)
     call sort_vec(scrtch,lms_rates_table(1)%nye,1,hydro_quantities(i:i)%ye,indx_yel(i:i),1,1)
     scrtch=lms_rates_table(1)%ye(:,idn1)
     call sort_vec(scrtch,lms_rates_table(1)%nye,1,hydro_quantities(i:i)%ye,indx_yer(i:i),1,1)
  enddo

! do the linear interpolation in the table
  do i=1,nnmax

!        the following is done in order to get reasonable values
!        also beyond the extent of the table
     idn =max(indx_dn(i) ,1) ; idn1 =min(idn+1,lms_rates_table(1)%nro)
     itm =max(indx_tm(i) ,1) ; itm1 =min(itm+1,lms_rates_table(1)%ntt)
     iyel=max(indx_yel(i),1) ; iyel1=min(iyel+1,lms_rates_table(1)%nye)
     iyer=max(indx_yer(i),1) ; iyer1=min(iyer+1,lms_rates_table(1)%nye)


     if (idn1.gt.idn) then
        x_dn=(hydro_quantities(i)%rho-lms_rates_table(1)%rho(idn))/(lms_rates_table(1)%rho(idn1)-lms_rates_table(1)%rho(idn))

        if (config%lms_spectra) then
           x_ld=(log_den(i)-lms_rates_table(1)%log_rho(idn))/(lms_rates_table(1)%log_rho(idn1)-lms_rates_table(1)%log_rho(idn))
        endif

     else
        x_dn=0._rk

        if (config%lms_spectra) x_ld=0._rk
     endif


     if (itm1.gt.itm) then
        x_tm=(hydro_quantities(i)%tem-lms_rates_table(1)%tem(itm))/(lms_rates_table(1)%tem(itm1)-lms_rates_table(1)%tem(itm))
     else
        x_tm=0._rk
     endif

     if (iyel1.gt.iyel) then
        x_yel=(hydro_quantities(i)%ye-lms_rates_table(1)%ye(iyel,idn))/ &
              (lms_rates_table(1)%ye(iyel1,idn)-lms_rates_table(1)%ye(iyel,idn))
     else
        x_yel=0._rk
     endif
     if (iyer1.gt.iyer) then
        x_yer=(hydro_quantities(i)%ye-lms_rates_table(1)%ye(iyer,idn1))/ &
              (lms_rates_table(1)%ye(iyer1,idn1)-lms_rates_table(1)%ye(iyer,idn1))
     else
        x_yer=0._rk
     endif

!        compute a continuous measure for the "distance" to table
     if (in_lmstb(i) .eq. -1._rk) then
        dstd=0._rk ; dstt=0._rk ; dsty=0._rk
        if (hydro_quantities(i)%rho .lt. lms_rates_table(1)%romin .or. hydro_quantities(i)%rho .gt. lms_rates_table(1)%romax) then
           dstd=((hydro_quantities(i)%rho-lms_rates_table(1)%rho(idn))/den_bnd)**2
        endif
        if (hydro_quantities(i)%tem .lt. lms_rates_table(1)%ttmin .or. hydro_quantities(i)%tem .gt. lms_rates_table(1)%ttmax) then
           dstt=((hydro_quantities(i)%tem-lms_rates_table(1)%tem(itm))/tem_bnd)**2
        endif
        if (hydro_quantities(i)%ye .lt. lms_rates_table(1)%yemin .or. hydro_quantities(i)%ye .gt. lms_rates_table(1)%yemax) then
!          if (iyel.ne.iyer) write(*,*) 'PANIC:',i,iyel,iyer
           dsty=((hydro_quantities(i)%ye -lms_rates_table(1)%ye(iyel,idn))/yee_bnd)**2
        endif
        dist(i)=min(dstd+dstt+dsty,1._rk)**2
     else
        dist(i)=0._rk
     endif


! trilinear interpolation
#define ARRAY_SIZE 3
#define RATE rtot(i)
#define FIELD lms_rates_table(1)%rtot
#include "electron_capture_on_nuclei.interpolate_rates.X90"

     if (config%lms_spectra) then
        ! neutrino equilbrium chemical potential
#define ARRAY_SIZE 3
#define RATE chenu(i)
#define FIELD lms_rates_table(1)%mu_nu
#include "electron_capture_on_nuclei.interpolate_rates.X90"

        ! spectra (engery slice)
#define ARRAY_SIZE 4
#define RATE spectra(i,:)
#define FIELD lms_rates_table(1)%spectra
#include "electron_capture_on_nuclei.interpolate_rates.X90"

     endif

     if ( (config%lms_spectra .and. config%qfit) .or. .not.config%lms_spectra)  then
#define ARRAY_SIZE 3
#define RATE chel(i)
#define FIELD lms_rates_table(1)%chel
#include "electron_capture_on_nuclei.interpolate_rates.X90"

#define ARRAY_SIZE 3
#define RATE qeff(i)
#define FIELD lms_rates_table(1)%qeff
#include "electron_capture_on_nuclei.interpolate_rates.X90"
     endif

     if (.not.config%lms_spectra) then
#define ARRAY_SIZE 3
#define RATE yhvy(i)
#define FIELD lms_rates_table(1)%yhvy
#include "electron_capture_on_nuclei.interpolate_rates.X90"

        if (.not.config%lms_old) then
! new lms used
#define ARRAY_SIZE 3
#define RATE chep(i)
#define FIELD lms_rates_table(1)%chep
#include "electron_capture_on_nuclei.interpolate_rates.X90"

#define ARRAY_SIZE 3
#define RATE chen(i)
#define FIELD lms_rates_table(1)%chen
#include "electron_capture_on_nuclei.interpolate_rates.X90"
        endif ! lms_old
     endif ! .not. spectra

  enddo

  if (config%lms_spectra) then
     rtot=exp(rtot)
     if (.not.config%qfit) spectra=exp(spectra)
  endif

  do i=1,nnmax
! FFN-scratch
     fnp=max(hydro_quantities(i)%zNuc-20._rk , 0._rk)
     fnp=min(fnp , 8._rk)
     fnh=max(40._rk-(hydro_quantities(i)%aNuc-hydro_quantities(i)%zNuc) , 0._rk)
     fnh=min(fnh , 6._rk)
     rffn(i)=hydro_quantities(i)%q_ga**2 *fnp*fnh*wc_s0/(14._rk*wc_me**2)*wc_ga2*hydro_quantities(i)%dnuc

! LMS-scratch
! beware rtot must be multiplied by the heavy NUMBER fraction
! this is done here

     if (.not.config%lms_old) then

        if (config%lms_heavy) then
           ! use number fraction of LMS-table (folding was done by initialisation)
           rtot(i)=hydro_quantities(i)%q_ga**2 *rtot(i)*cnst_phsp*(hydro_quantities(i)%rho/pc_mb)
        else

! use the NUMBER fraction of EOS:
           if (hydro_quantities(i)%dnuc .eq. 0._rk) then
              rtot(i)=0._rk
           else
              rtot(i)=hydro_quantities(i)%q_ga**2 *rtot(i)*cnst_phsp*hydro_quantities(i)%dnuc
           endif

        endif ! lms_heavy

     else ! old_lms
        rtot(i)=hydro_quantities(i)%q_ga**2 * rtot(i)*cnst_phsp*(hydro_quantities(i)%rho/pc_mb)
     endif

  enddo

! energy dependence proposed by LMS is the same as FFN,
!     only q-value is different

        do ie=1,nemax

           do i=1,nnmax


              if (.not.config%lms_spectra .or. (config%lms_spectra .and. &
                                                config%qfit) ) then
                 cnu_m = hydro_quantities(i)%celectron+hydro_quantities(i)%cproton-hydro_quantities(i)%cneutron - wc_mq
              !            cnu_m = chel(i)+chep(i)-chen(i) - wc_mq
                 eargn = beta(i)*(em(ie)      -cnu_m)

                 ! LMS:
                 ! if qeff > -0.511 MeV (electron restmass) than the neutrino-energy
                 ! will be larger than the initial electron energy
                 ! However, since the electron energy must be atleast the restmass
                 ! this imposes the constraint that the minium neutrino energy
                 ! must be max(restmass +q,0) for qeff .gt. -0.511 MeV
                 ! if this is not the case the spectra must be zero.

                 earge = beta(i)*(em(ie)-qeff(i)-chel(i))

                 if (.not.config%lms_old) then

                    if ((qeff(i) .lt. -wc_me) .or.             &
                         (qeff(i) .gt. -wc_me .and. em(ie) .ge. &
                         wc_me+qeff(i))) then

                       if (eargn.gt.35._rk .and. earge.gt.35._rk) then
                          th = exp(min(eargn-earge,bigeexp))
                       else
                          th = ( 1._rk + exp(min(eargn,bigeexp)) ) / &
                               ( 1._rk + exp(min(earge,bigeexp)) )
                       endif

                       if (hydro_quantities(i)%dnuc .eq. 0._rk) then
                          rlms_e=0._rk
                       else

                          rlms_e = (em(ie)-qeff(i))**2*th*rtot(i)
                       endif

                    else
                       rlms_e=0._rk
                    endif
                 else !.not.config_lms_old

                       if (eargn.gt.35._rk .and. earge.gt.35._rk) then
                          th = exp(min(eargn-earge,bigeexp))
                       else
                          th = ( 1._rk + exp(min(eargn,bigeexp)) ) / &
                               ( 1._rk + exp(min(earge,bigeexp)) )
                       endif

                       if (hydro_quantities(i)%dnuc .eq. 0._rk) then
                          rlms_e=0._rk
                       else

                          rlms_e = (em(ie)-qeff(i))**2*th*rtot(i)
                       endif

                 endif !  .not.config_lms_old

           else ! .not.config%lms_spectra

              ! interpolate spectra

              !c LMS:
              ! we still need the neutrino equilibrium distribution
              ! (\kappa=j\f_eq)

              earge = beta(i)*(em(ie) - chenu(i))
              cnu_m = hydro_quantities(i)%celectron+hydro_quantities(i)%cproton-hydro_quantities(i)%cneutron - wc_mq
              eargn = beta(i)*(em(ie) - cnu_m)

          if (eargn.gt.35._rk.and.earge.gt.35._rk) then
             th = exp(min(eargn-earge,bigeexp))
          else
             th = (1.0_rk + exp(min(eargn,bigeexp)))/ &
                  (1.0_rk + exp(min(earge,bigeexp)))
          endif

          if (hydro_quantities(i)%dnuc .eq. 0._rk) then
             rlms_e=0._rk
          else

!     j=n(\epsilon)/(\epsilon)^2
             rlms_e = spectra (i,ie) * th * rtot (i) / em (ie) ** 2
          endif


       endif ! .not.config%lms_spectra
! FFN:
              cnu_m = hydro_quantities(i)%celectron+hydro_quantities(i)%cproton-hydro_quantities(i)%cneutron - wc_mq ! new
              eargn= beta(i)*(em(ie)      -cnu_m)  ! new

              qffn=-((hydro_quantities(i)%cneutron-hydro_quantities(i)%cproton+wc_mq) + delta)
              earge = beta(i)*(em(ie)-qffn-hydro_quantities(i)%celectron)


              if (eargn.gt.35._rk .and. earge.gt.35._rk) then
                 th = exp(min(eargn-earge,bigeexp))
              else
                 th = ( 1._rk + exp(min(eargn,bigeexp)) ) / &
                      ( 1._rk + exp(min(earge,bigeexp)) )
              endif

              rffn_e = (em(ie)-qffn)**2*th*rffn(i)

! weighted interpolation between LMS-rates and FFN-rates,
!   for this purpose we assume m_e<<e_neutrino
!              if (i .eq. 1. .and. ie .eq. 2) then
!                 print *,"new ",opacity(i,ie)%aab, dist(i),rlms_e,rffn_e
!              endif

              aab(i,ie) = (1._rk-dist(i))*rlms_e+dist(i)*rffn_e


!              if (i .eq. 1. .and. ie .eq. 2) then
!                 print *,"new ",opacity(i,ie)%aab, dist(i),rlms_e,rffn_e
!              endif


           enddo
        enddo


!  call timer%stop("abso_nuci_lms")
  return

end function abso_nuci_lms


end module lms_rates

module abso_nuci_ffn_rates

  implicit none
  private

  public abso_nuci_ffn


  contains
! -------------------------------------------------------------------

    function abso_nuci_ffn(em, isi, isf, hydro_quantities) result(aab)


! Purpose: calculate opacities for
!            absorption of neutrinos on nuclei
!
! Input:
!        emid:        Neutrino energies
!        aab :        absorption opacity (due to other interactions)
!        dnuc:        NUMBER density of nuclei AZ,AN,AA
!        az,aa:       Z,A
!        cel :        chem potential of electrons (including restmass)
!        cne :        chem potential of neutrons  (without restmass)
!        cpr :        chem potential of protons   (without restmass)
!        aab :        opacity due to other processes
!
! Output:
!        aab:         aab + opacities for abs. of neutrinos on nuclei
!
! Reference: Bruenn ApJSS 58, 771 (1985) Appendix C
!
      use precision

      use phycon
      use machcons
      use abort
      use timings

      use mod_rate_types, only : hydro_quantities_t, absorption_opacity_t

      use configure
      implicit none
! LOCAL variables that are not in modules
      integer(kind=ik), intent(in)           :: isi, isf
      real(kind=rk), intent(in)              :: em(:)
      type(hydro_quantities_t), intent(in)   :: hydro_quantities(:)
      real(kind=rk), dimension(size(hydro_quantities), &
                size(em))                :: aab
!      type(absorption_opacity_t),       &
!      dimension(size(hydro_quantities), &
!                size(em))                    :: opacity


      integer(kind=ik)                       :: i, ie, nemax, nnmax
      real(kind=rk)                          :: con_a, th, earge, eargn, cel_m, &
                                                cnu_m, qprim
      real(kind=rk)                          :: an, fnh, fnp

      real(kind=rk), &
      dimension(size(hydro_quantities)) :: beta, fn
      real(kind=rk) ,parameter         :: delta=3._rk    !Excitation energy [MeV]


!      call timer%start("abso_nuci_ffn")

      nnmax=size(hydro_quantities)
      nemax=size(em)

      aab(:,:) = 0.0_rk

      do i=1,nnmax
         beta(i)=1._rk/hydro_quantities(i)%tem

!         if (i .eq. 1) then
!            print *,"new a",hydro_quantities(i)%zNuc,hydro_quantities(i)%aNuc
!         endif


         if (hydro_quantities(i)%zNuc .le. 20._rk) then
            fnp=0.0
         elseif (hydro_quantities(i)%zNuc .gt. 20._rk .and. hydro_quantities(i)%zNuc .le. 28._rk) then
            fnp=hydro_quantities(i)%zNuc-20._rk
         else
            fnp=8._rk
         endif

         an=hydro_quantities(i)%aNuc-hydro_quantities(i)%zNuc
         if (an .le. 34._rk) then
            fnh=6._rk
         elseif(an .gt. 34._rk .and. an .le. 40._rk) then
            fnh=40._rk-an
         else
            fnh=0._rk
         endif

         fn(i)=fnp*fnh
      enddo

      do ie=1,nemax

         do i=1,nnmax
            qprim=(hydro_quantities(i)%cneutron -hydro_quantities(i)%cproton+wc_mq) + delta

            cnu_m = hydro_quantities(i)%celectron+hydro_quantities(i)%cproton-hydro_quantities(i)%cneutron - wc_mq
            cel_m = hydro_quantities(i)%celectron

            eargn = beta(i)*(em(ie)      -cnu_m)
            earge = beta(i)*(em(ie)+qprim-cel_m)

            if (eargn.gt.35._rk .and. earge.gt.35._rk) then
               th = exp(min(eargn-earge,bigeexp))
            else
               th = ( 1._rk + exp(min(eargn,bigeexp)) ) / &
                    ( 1._rk + exp(min(earge,bigeexp)) )
            endif

            if (em(ie)+qprim .le. wc_me) then
               con_a = 0._rk
            else
               con_a = 1._rk/14._rk*wc_s0*wc_ga2*hydro_quantities(i)%q_ga**2 &
                    * (em(ie)+qprim)/wc_me**2 &
                    * sqrt((em(ie)+qprim)**2-wc_me**2)
            endif
!            if (i .eq. 1 .and. ie .eq. 1 ) then
!               print *,"new ",i,ie,opacity(i,ie)%aab
!            endif
            aab(i,ie) = con_a*th*hydro_quantities(i)%dnuc*fn(i)
        enddo
      enddo

!      call timer%stop("abso_nuci_ffn")
      return

    end function abso_nuci_ffn

end module abso_nuci_ffn_rates
