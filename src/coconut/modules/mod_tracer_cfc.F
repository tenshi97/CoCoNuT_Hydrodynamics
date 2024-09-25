c     ============================================================
      MODULE tracer_cfc
c     ============================================================
      USE precision
      USE size_cfc, ONLY: m, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc
      use configure

      IMPLICIT NONE

      SAVE

      integer (kind=ik) :: tracer_num
      real (kind=rk), allocatable :: pos_tracer_r(:)
      real (kind=rk), allocatable :: pos_tracer_r_old(:)
      real (kind=rk), allocatable :: pos_tracer_theta(:)
      real (kind=rk), allocatable :: pos_tracer_theta_old(:)
      real (kind=rk), allocatable :: v_1_tracer(:)
      real (kind=rk), allocatable :: v_2_tracer(:)
      real (kind=rk), allocatable :: temp_tracer(:)
      real (kind=rk), allocatable :: xnnu_tracer(:,:)
      real (kind=rk), allocatable :: den_tracer(:)
      real (kind=rk), allocatable :: pressure_tracer(:)
      real (kind=rk), allocatable :: entropy_tracer(:)
      real (kind=rk) :: tracer_printtime = -1
      integer (kind=ik), allocatable :: tracer_ci_r(:)
      integer (kind=ik), allocatable :: tracer_ci_theta(:)
      LOGICAL, allocatable :: neg_check(:)
      LOGICAL, allocatable :: nan_check(:)
      integer (kind=ik), allocatable :: tracer_id(:)
      integer (kind=ik), allocatable ::tracer_proc(:)
      CONTAINS
c     ============================================================

      SUBROUTINE allocate_tracer_cfc
        
c     ============================================================
      USE mo_mpi
      USE precision
      USE size_cfc on
      IMPLICIT NONE
      
      integer :: istat
      integer(kind=ik) :: i
      character*64 :: dfilename,procstr,outfilename,mpifilename
      character*64 :: suffix
      suffix = '.dat'
      WRITE(unit=procstr,fmt='(i4)') myproc
      procstr = trim(adjustl(procstr))
      outfilename = 'debugmpi/filen'//procstr // suffix
      dfilename = 'debugmpi/debugn'// procstr // suffix
      mpifilename = 'debugmpi/mpin' // procstr // suffix
c      WRITE(unit=dfilename,fmt="('debugmpi/file',i,'.dat')") myproc
      OPEN(UNIT=7777,FILE=outfilename)
      OPEN(UNIT=8888,FILE=dfilename)
      open(UNIT=9999,FILE=mpifilename)
      tracer_num = 20
      allocate(
     &  pos_tracer_r(tracer_num),
     &  pos_tracer_r_old(tracer_num),
     &  pos_tracer_theta(tracer_num),
     &  pos_tracer_theta_old(tracer_num),
     &  v_1_tracer(tracer_num),
     &  v_2_tracer(tracer_num),
     &  temp_tracer(tracer_num),
     &  xnnu_tracer(tracer_num,1:config%qn),
     &  den_tracer(tracer_num),
     &  pressure_tracer(tracer_num),
     &  entropy_tracer(tracer_num), 
     &  tracer_ci_r(tracer_num),
     &  tracer_ci_theta(tracer_num),
     &  nan_check(tracer_num),
     &  neg_check(tracer_num),
     &  tracer_id(tracer_num),
     &  tracer_proc(tracer_num),
     &  stat=istat)
      tracer_ci_r = -2_ik
      tracer_ci_theta =-2_ik
      neg_check = .true.
      nan_check = .true.
      do i=1,tracer_num
            tracer_id(i) = myproc*tracer_num+i
      end do
      return
c     xnnu 1=neutron 2=proton 3= D/ 4= T  5=3_He 6= 4_He....
      END SUBROUTINE allocate_tracer_cfc
c     ==================================================================

      END MODULE tracer_cfc

c     ==================================================================
