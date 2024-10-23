c     ============================================================
      MODULE tracer_cfc
c     ============================================================
      USE precision
      USE size_cfc, ONLY: m, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc
      use configure

      IMPLICIT NONE

      SAVE

      integer (kind=ik) :: tracer_num,tracer_pp,tracer_total
      real (kind=rk), allocatable :: pos_tracer_r(:)
      real (kind=rk), allocatable :: pos_tracer_r_old(:)
      real (kind=rk), allocatable :: pos_tracer_theta(:)
      real (kind=rk), allocatable :: pos_tracer_theta_old(:)
      real (kind=rk), allocatable :: v_1_tracer(:)
      real (kind=rk), allocatable :: v_2_tracer(:)
      real (kind=rk) :: tracer_printtime = -1
      integer (kind=ik), allocatable :: tracer_ci_r(:)
      integer (kind=ik), allocatable :: tracer_ci_theta(:)
      integer (kind=ik), allocatable :: tracer_id(:)
      integer (kind=ik), allocatable ::tracer_proc(:)
      integer (kind=ik) :: qyrays_t
      integer (kind=ik) n_ts,n_te
      CONTAINS
c     ============================================================

      SUBROUTINE allocate_tracer_arrays(tracer_rstrt)
        
c     ============================================================
      USE mo_mpi
      USE precision
      USE size_cfc
c      USE gr_input, only: trid_rstrt, trx_rstrt, try_rstrt, 
c     &  trnum_rstrt
      IMPLICIT NONE
      
      integer :: istat
      integer(kind=ik) :: i,j
      character*64 :: dfilename,procstr,outfilename,mpifilename
      character*64 :: suffix
      logical, intent(in) :: tracer_rstrt
      suffix = 'dat'
      WRITE(unit=procstr,fmt='(i4)') myproc
      procstr = trim(adjustl(procstr))
      outfilename = 'debugmpi/file'//procstr // '.'//suffix
      dfilename = 'debugmpi/debug'// procstr // '.'//suffix
      mpifilename = 'debugmpi/mpi' // procstr // '.'//suffix
      OPEN(UNIT=7777,FILE=outfilename)
      OPEN(UNIT=8888,FILE=dfilename)
      open(UNIT=9999,FILE=mpifilename)
      if(tracer_rstrt) then
        return 
      end if
      tracer_pp = config%trpp
      tracer_num = tracer_pp
      tracer_total = tracer_pp * nprocs
      qyrays_t = n_e-n_s+1
      WRITE(7777,*) "INIT INFO:tracer per process:",tracer_pp
      WRITE(7777,*) "INIT INFO:tracer total",tracer_total
      WRITE(7777,*) "INIT INFO:y-rays total",qyrays_t
      allocate(
     &  pos_tracer_r(tracer_num),
     &  pos_tracer_r_old(tracer_num),
     &  pos_tracer_theta(tracer_num),
     &  pos_tracer_theta_old(tracer_num),
     &  v_1_tracer(tracer_num),
     &  v_2_tracer(tracer_num),
     &  tracer_ci_r(tracer_num),
     &  tracer_ci_theta(tracer_num),
     &  tracer_id(tracer_num),
     &  tracer_proc(tracer_num),
     &  stat=istat)
      tracer_ci_r = -2_ik
      tracer_ci_theta =-2_ik
      n_ts = myproc*tracer_num+1
      n_te = myproc*tracer_num+tracer_num
      do i=1,tracer_num
        tracer_id(i) = myproc*tracer_num+i
      end do

      return
c     xnnu 1=neutron 2=proton 3= D/ 4= T  5=3_He 6= 4_He....
      END SUBROUTINE allocate_tracer_arrays
c     ==================================================================
      END MODULE tracer_cfc
