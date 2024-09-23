module data_lms_rates

  use precision
  use configure

  implicit none

  private
  public load_lmstable, iye_tb


  integer(kind=ik) :: idn
  integer(kind=ik) :: idn_tb=14, itm_tb=12, iye_tb=12, ien_tb=0

  real(kind=rk)    :: dmn_tb,dmx_tb,tmn_tb,tmx_tb,ymn_tb,ymx_tb


  real(kind=rk), allocatable, target ::  den_tb(:),tem_tb(:), &
                                         yee_tb(:,:), log_den_tb(:)


  real(kind=rk),  dimension(:,:,:), &
                   allocatable, target  :: chel_tb,yp_tb,chep_tb,yn_tb,&
                                           chen_tb,zav_tb,aav_tb,zh_tb,&
                                           ah_tb,yhvy_tb,yalpha_tb,    &
                                           rtot_tb,rnrm_tb,qeff_tb,    &
                                           mu_nu_tb,n_tb,mu_e_tb
  real(kind=rk),  dimension(:,:,:,:), &
                   allocatable, target  :: spectra_tb

  real(kind=rk), dimension (:), allocatable, target :: energy_tb

  character(LEN=30) ::  tablefile
contains

  subroutine load_lmstable
! LOCAL variables that are not in modules
    use precision
    use rates_data_type, only : lms_rates_table

    use energy_grid_rt

    use abort
    use mpi_vertex, only : myproc
    use print_stdout_mod

    use configure
    implicit none

    integer(kind=ik) :: idn_rd,itm_rd,iye_rd,ien_rd
    integer(kind=ik) :: istat


    integer(kind=ik) :: i,j,k
#ifdef CRAY
    logical          :: ierr
#endif

    if (config%lms_old) then
       tablefile="./tables/lmsscreened-nsetable"
    endif


    if (.not.config%lms_old .and. .not.config%lms_spectra) then
       idn_tb=14
       itm_tb=12
       iye_tb=12
       tablefile="./tables/lmsscreened-heavy"
    endif


    if (config%lms_spectra) then
       idn_tb=25
       itm_tb=35
       iye_tb=13
       ien_tb=config%iemax
        tablefile="./tables/spektrale_lms_tabelle"
    endif

#ifdef CRAY
    call asnunit(53,' -F f77 -N ieee',ierr)
#endif

      allocate(den_tb(idn_tb),tem_tb(itm_tb),yee_tb(iye_tb,idn_tb),stat=istat)

      if (istat /= 0) raise_abort("load_lmstable(): error in allocating arrays")

      allocate(chel_tb(idn_tb,itm_tb,iye_tb),                            &
           yp_tb(idn_tb,itm_tb,iye_tb),chep_tb(idn_tb,itm_tb,iye_tb),    &
           yn_tb(idn_tb,itm_tb,iye_tb),chen_tb(idn_tb,itm_tb,iye_tb),    &
           zav_tb(idn_tb,itm_tb,iye_tb),aav_tb(idn_tb,itm_tb,iye_tb),    &
           zh_tb(idn_tb,itm_tb,iye_tb),ah_tb(idn_tb,itm_tb,iye_tb),      &
           yhvy_tb(idn_tb,itm_tb,iye_tb),yalpha_tb(idn_tb,itm_tb,iye_tb),&
           rtot_tb(idn_tb,itm_tb,iye_tb),rnrm_tb(idn_tb,itm_tb,iye_tb),  &
           qeff_tb(idn_tb,itm_tb,iye_tb),stat=istat)
      if (istat /= 0) raise_abort("load_lmstable(): error in allocating arrays")

      if (config%lms_spectra) then
         allocate(energy_tb(ien_tb),spectra_tb(idn_tb,itm_tb,iye_tb,ien_tb), &
                mu_nu_tb(idn_tb,itm_tb,iye_tb),                            &
                mu_e_tb(idn_tb,itm_tb,iye_tb), stat=istat)

         if (istat /= 0) raise_abort("load_lmstable(): error in allocating arrays")

      endif

      open(53,file=tablefile,form = 'unformatted')

      if (.not.config%lms_spectra) then
         read(53) idn_rd,itm_rd,iye_rd
      else
         read(53) idn_rd,itm_rd,iye_rd,ien_rd
      endif

      if(idn_rd.ne.idn_tb.or.itm_rd.ne.itm_tb.or.iye_rd.ne.iye_tb) then
         write(*,'(" idn_tb = ",i4," idn_rd = ",i4)') idn_tb,idn_rd
         write(*,'(" itm_tb = ",i4," itm_rd = ",i4)') itm_tb,itm_rd
         write(*,'(" iye_tb = ",i4," iye_rd = ",i4)') iye_tb,iye_rd
         if (config%lms_spectra) then
            write(*,'(" ien_tb = ",i4," ien_rd = ",i4)') ien_tb,ien_rd
         endif
         raise_abort("load_lmstable(): incorrect dimensions of table")
      endif

      if (.not.config%lms_old .and. .not.config%lms_spectra) then
         read(53) den_tb,tem_tb,yee_tb,chel_tb,yp_tb,chep_tb,yn_tb, &
               chen_tb,zav_tb,aav_tb,zh_tb,ah_tb,yhvy_tb,        &
               yalpha_tb,rtot_tb,qeff_tb,rnrm_tb
      endif


      if (config%lms_old) then
         
         read(53) den_tb,tem_tb,yee_tb
         read(53) qeff_tb,chel_tb,rtot_tb,yhvy_tb,rnrm_tb
      endif

      if (config%lms_spectra) then

         read(53) den_tb,tem_tb,yee_tb,energy_tb,yhvy_tb,mu_nu_tb
         if (maxval(abs(emid(1:config%iemax)-energy_tb(1:config%iemax))) .gt.1.0e-10_rk) then
            call printit_taskX(0,"emid:" , emid)
            call printit_taskX(0,"energy_tb: ", energy_tb)
            raise_abort("load_lmstable(): energy grids incompatible")
         else
            call printit_taskX(0,"Energy grid correct.")
         end if

         read(53) spectra_tb, rtot_tb

         if (config%qfit) read(53) n_tb, qeff_tb, mu_e_tb

      log_den_tb = log (den_tb)
   endif


!      read(53) den_tb,tem_tb,yee_tb
!      read(53) qeff_tb,chel_tb,rtot_tb,yhvy_tb,rnrm_tb
       close(53)

       if (myproc .eq. 0) then
          write (*,*) idn_tb,itm_tb,iye_tb
       endif

       if (.not.config%lms_spectra) then
          ! normalisation but still multiplication with yhvy necessary

          if (myproc .eq. 0) then
             write (*,*) size(rtot_tb,dim=1),size(rtot_tb,dim=2),size(rtot_tb,dim=3)
             write (*,*) size(rnrm_tb,dim=1),size(rnrm_tb,dim=2),size(rnrm_tb,dim=3)
          endif

          do i=1,14
             do j=1,12
                do k=1,12
                   if (rnrm_tb(i,j,k) .eq. 0._rk) then
                      write (*,*) i,j,k,rnrm_tb(i,j,k)
                   endif
                enddo
             enddo
          enddo

          where( rnrm_tb .ne. 0._rk)
             rtot_tb=rtot_tb/rnrm_tb
          endwhere
          

          if (.not.config%lms_old .and. config%lms_heavy) then
             ! we do already here the folding with Y_heavy of the LMS-Table
             rtot_tb = rtot_tb * yhvy_tb
          endif
       else ! config%lms_spectra


          if (config%lms_heavy) then
             ! we do already here the folding with Y_heavy of the LMS-Table
             rtot_tb = rtot_tb * yhvy_tb
          endif

           if (config%qfit) then
              ! we also need the normalization factor for the q-fit
              rtot_tb = rtot_tb * n_tb
           endif


           rtot_tb = log (rtot_tb)
           spectra_tb = log (spectra_tb)
        endif ! config%lms_spectra



        if (.not.config%lms_old .and. .not.config%lms_spectra) then
           deallocate(yp_tb,yn_tb,zav_tb,aav_tb,zh_tb,ah_tb,yalpha_tb)
        endif

        if (config%lms_old) then
           deallocate(yp_tb,yn_tb,zav_tb,aav_tb,zh_tb,ah_tb,yalpha_tb, &
                  chep_tb,chen_tb)
        endif

!       chen_tb(:,:,:)=1.
!       chep_tb(:,:,:)=1.
!       yhvy_tb(:,:,:)=1.



       dmn_tb=den_tb(1)   ; dmx_tb=den_tb(idn_tb)
       tmn_tb=tem_tb(1)   ; tmx_tb=tem_tb(itm_tb)
       ymn_tb=yee_tb(1,1) ; ymx_tb=yee_tb(iye_tb,1)


! assign pointers to data_type

       lms_rates_table(1)%name          = tablefile

       if (config%lms_spectra) then
          lms_rates_table(1)%nene           = ien_tb
       else
          lms_rates_table(1)%nene           = 0
       endif

       lms_rates_table(1)%nro           = idn_rd
       lms_rates_table(1)%ntt           = itm_rd
       lms_rates_table(1)%nye           = iye_rd

       lms_rates_table(1)%romax         = dmx_tb
       lms_rates_table(1)%romin         = dmn_tb
       lms_rates_table(1)%ttmax         = tmx_tb
       lms_rates_table(1)%ttmin         = tmn_tb
       lms_rates_table(1)%yemax         = ymx_tb
       lms_rates_table(1)%yemin         = ymn_tb


       lms_rates_table(1)%rho           => den_tb
       lms_rates_table(1)%tem           => tem_tb
       lms_rates_table(1)%ye            => yee_tb

       if (config%lms_spectra) then
          lms_rates_table(1)%n             => n_tb
          lms_rates_table(1)%log_rho       => log_den_tb
          lms_rates_table(1)%spectra       => spectra_tb
          lms_rates_table(1)%energy        => energy_tb
          lms_rates_table(1)%mu_nu         => mu_nu_tb
          lms_rates_table(1)%mu_e          => mu_e_tb
       else
          
          lms_rates_table(1)%chel          => chel_tb
          lms_rates_table(1)%chep          => chep_tb
          lms_rates_table(1)%chen          => chen_tb
       endif

       lms_rates_table(1)%yhvy          => yhvy_tb
       lms_rates_table(1)%rtot          => rtot_tb
       lms_rates_table(1)%rnrm          => rnrm_tb
       lms_rates_table(1)%qeff          => qeff_tb

       if (.not.config%lms_old .and. .not.config%lms_spectra) then
          nullify(lms_rates_table(1)%yp,lms_rates_table(1)%yn,   &
               lms_rates_table(1)%zav,lms_rates_table(1)%aav, &
               lms_rates_table(1)%zh,lms_rates_table(1)%ah,   &
               lms_rates_table(1)%yalpha,                     &
               lms_rates_table(1)%log_rho,                &
               lms_rates_table(1)%mu_nu,lms_rates_table(1)%mu_e, &
               lms_rates_table(1)% energy)
       endif

       if (config%lms_old) then
             nullify(lms_rates_table(1)%yp,lms_rates_table(1)%yn,   &
               lms_rates_table(1)%zav,lms_rates_table(1)%aav,       &
               lms_rates_table(1)%zh,lms_rates_table(1)%ah,         &
               lms_rates_table(1)%yalpha,lms_rates_table(1)%chep,   &
               lms_rates_table(1)%chen,                             &
               lms_rates_table(1)%log_rho,                      &
               lms_rates_table(1)%mu_nu,lms_rates_table(1)%mu_e, &
               lms_rates_table(1)% energy)
          endif

          if (config%lms_spectra) then
             nullify(lms_rates_table(1)%yp,lms_rates_table(1)%yn,        &
               lms_rates_table(1)%zav,lms_rates_table(1)%aav,     &
               lms_rates_table(1)%zh,lms_rates_table(1)%ah,       &
               lms_rates_table(1)%yalpha,lms_rates_table(1)%chel, &
               lms_rates_table(1)%chep,lms_rates_table(1)%chen)

          endif


       call printit_taskX(0,"================================================================================")
       call printit_taskX(0," ")
       call printit_taskX(0,"   load_lmstable> table --> ",tablefile)
       call printit_taskX(0,"   <-- successfully installed!")
       call printit_taskX(0,"   idn_tb = ",lms_rates_table(1)%nro)
       call printit_taskX(0,"   itm_tb = ",lms_rates_table(1)%ntt)
       call printit_taskX(0,"   iye_tb = ",lms_rates_table(1)%nye)
       call printit_taskX(0,"   romin  = ",lms_rates_table(1)%romin)
       call printit_taskX(0,"   romax  = ",lms_rates_table(1)%romax)
       call printit_taskX(0,"   ttmin  = ",lms_rates_table(1)%ttmin)
       call printit_taskX(0,"   ttmax  = ",lms_rates_table(1)%ttmax)
       call printit_taskX(0,"   yemin  = ",lms_rates_table(1)%yemin)
       call printit_taskX(0,"   yemax  = ",lms_rates_table(1)%yemax)
       call printit_taskX(0," ")
       call printit_taskX(0,"   WARNING: *** the table provided by LMS of was extrapolated ***")
       call printit_taskX(0,"   original boundaries of ye(rho) are")
       call printit_taskX(0,"   irho  rho [g/ccm] yemin yemax")

       do idn=1,idn_tb
          if (myproc .eq. 0) then
             write(*,'(3x,I4,1pe12.3,2e12.3)') idn, &
                       & lms_rates_table(1)%rho(idn),lms_rates_table(1)%ye(2,idn),lms_rates_table(1)%ye(iye_tb-1,idn)
          endif
       enddo

       call printit_taskX(0," ")
       if (.not.config%lms_old) then
          if (.not.config%lms_heavy) then
             call printit_taskX(0,"lmsrates: Using Y_heavy from high density EOS table")
          else
             call printit_taskX(0,"lmsrates: Using Y_heavy from LMS table")
          endif
       endif

      if (config%lms_old) then
         call printit_taskX(0,"lmsrates: Y_heavy of LMS table is used (lms_old is set)")
      endif

       call printit_taskX(0,"================================================================================")


  end subroutine load_lmstable

end module data_lms_rates
