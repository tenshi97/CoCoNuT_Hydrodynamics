
c     ******************************************************************

      MODULE tracer

c     ******************************************************************

      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE configure
      IMPLICIT NONE

!     r contains cell-centred coordinats (grid runs from index 1 to index m)
!     r_if contains cell boundries
!     delta_r contains r_if_outer - r_if_inner
      
!     v_1 is the radial velocity component




      CONTAINS
      
c     ==================================================================

      SUBROUTINE tci_search (tracer_index)
      
c     ==================================================================
      USE precision
      USE grid_cfc
      UsE hydro_primitives_cfc

      USE tracer_cfc
c     tracer_ci_r(i)---> the i-th tracer's cell index initialized as -2
      INTEGER (kind=ik) tci_r
      INTEGER (kind=ik), INTENT(IN) :: tracer_index
      REAL (kind=rk) :: pos_r
      pos_r = pos_tracer_r(tracer_index)
      tci_r = tracer_ci_r(tracer_index)
      if(tci_r.lt.0) then
            tracer_ci_r(tracer_index) = tci_search_bisection(pos_r)
      else
            if(tci_r.gt.1) then
                  if(pos_r.lt.r_if(tci_r-1)) then
                        if(pos_r.gt.r_if(tci_r-2)) then
                              tracer_ci_r(tracer_index) = tci_r-1
                        else
                              tracer_ci_r(tracer_index) = 
     &                        tci_search_bisection(pos_r)
                              print *,"[Warning] Unexpected Over-cross 
     &                        for particle",tracer_index
                        end if
                        return
                  end if
            end if
            if(tci_r.lt.m) then
                  if(pos_r.gt.r_if(tci_r)) then
                        if(pos_r.lt.r_if(tci_r+1)) then
                              tracer_ci_r(tracer_index) = tci_r+1
                        else
                              tracer_ci_r(tracer_index) = 
     &                        tci_search_bisection(pos_r)
                              print *,"[Warning] Unexpected Over-cross 
     &                        for particle",tracer_index 
                        end if
                        return
                  end if
            end if 
      end if
      END SUBROUTINE tci_search


c     ==================================================================

      FUNCTION tci_search_bisection (r_tracer)
      
c     ==================================================================

      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc

      USE tracer_cfc

      IMPLICIT NONE
      
      INTEGER (KIND=ik) tci_search_bisection
      REAL (KIND=rk), INTENT(IN) :: r_tracer
      INTEGER (KIND=ik) :: startPointer, endPointer, midPointer
      startPointer = 1
      endPointer = m
      IF(r_tracer<=r_if(0).or.r_tracer>=r_if(m)) THEN
            tci_search_bisection = -1
            RETURN
      END IF
      FINDER_LOOP: DO WHILE (startPointer<endPointer)
                        midPointer = (startPointer + endPointer) / 2
                        IF(r_tracer<r_if(midPointer)) THEN
                              endPointer = midPointer
                        ELSE
                              startPointer = midPointer+1
                        END IF
                  END DO FINDER_LOOP
                  tci_search_bisection = startPointer
                  
      END FUNCTION tci_search_bisection


c     ==================================================================

      SUBROUTINE tci_t_search (ti)
      
c     ==================================================================

      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc

      USE tracer_cfc

      IMPLICIT NONE
      
      INTEGER (KIND=ik), INTENT(IN) :: ti
      REAL (KIND=rk) :: theta_tracer
      INTEGER (KIND=ik) :: startPointer, endPointer, midPointer
      startPointer = n_s
      endPointer = n_e
      theta_tracer =pos_tracer_theta(ti)
      IF(theta_tracer<=theta_if(n_s-1)) THEN
         tracer_ci_theta(ti) = n_s-1
         RETURN
      END IF
      IF(theta_tracer>=theta_if(n_e)) THEN
         tracer_ci_theta(ti) = n_e+1
         RETURN
      END IF
      FINDER_LOOP: DO WHILE (startPointer<endPointer)
                        midPointer = (startPointer + endPointer) / 2
                        IF(theta_tracer<theta_if(midPointer)) THEN
                           endPointer = midPointer
                        ELSE
                           startPointer = midPointer+1
                        END IF
                  END DO FINDER_LOOP
                  tracer_ci_theta(ti) = startPointer

      END SUBROUTINE tci_t_search




c     ==================================================================

      SUBROUTINE tracer_position_init()
      
c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE size_cfc, ONLY: m, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc
      USE mo_mpi
      IMPLICIT NONE
      integer (KIND=ik) :: i,t_c
      t_c = n_s
      do i =1,m
         WRITE(7777,*) 'Radial Cell',i,'(',r_if(i-1),',',
     &   r_if(i),')'
      end do
      do i = n_s,n_e
            WRITE(7777,*) 'Theta Cell',i, '(', theta_if(i-1),',',
     &     theta_if(i),')'
      end do
      init_loop: do i=1,tracer_num
                  pos_tracer_r(i) = 1.5_rk**(i)*1e6;
                  pos_tracer_theta(i) = theta(t_c)
                  tracer_proc(i) = myproc
                  t_c = t_c + 1
                  if(t_c.gt.n_e) then
                        t_c = n_s
                  end if
                 call tci_search(i)
                 call tci_t_search(i)
      end do init_loop
      END SUBROUTINE tracer_position_init
c     Let's just start from 3rd-Runge-Kutta Method to predict the position evolution of the tracers, but first we should get the covariant velocity


c     ==================================================================

      SUBROUTINE tracer_integral(rkint_order,alpha_rk,beta_rk,delta_t)

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE size_cfc, ONLY: m, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc
      USE mo_mpi
      USE gfloat_hy, ONLY: time
      use intgrs_hy, only : nstep
      IMPLICIT NONE
      
      integer (KIND=ik) :: i,ir,it
      real (KIND=rk) :: v_1_t,v_2_t,temp_t      
      integer (KIND=ik), INTENT(IN) :: rkint_order
      real(KIND=rk), INTENT(IN) :: alpha_rk,beta_rk,delta_t
      WRITE(9999,*) "INTEGRATION START at time:",time," ",nstep
      LOOP_INT:DO i=1,tracer_num
            call tci_search(i)
            call tci_t_search(i)
            ir = tracer_ci_r(i)
            it = tracer_ci_theta(i)
            IF((ir.lt.0).or.(ir.gt.m))   THEN
                  CYCLE LOOP_INT
            END IF
            v_1_t = bilineary_interpolation
     &      (v_1(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i,
     &      .false.,"v_1")
     
            v_2_t = bilineary_interpolation
     &      (v_2(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i,
     &      .false.,"v_2")
            IF((it.lt.n_s))      THEN
               it = n_s
            ELSE IF(it.gt.n_e)   THEN
               it = n_e
            END IF

            v_1_t = g_up_11(ir,it,1)*v_1_t - 
     &      beta_up_1(ir,it,1)/alpha(ir,it,1)
            v_2_t = g_up_22(ir,it,1)*v_2_t-
     &      beta_up_2(ir,it,1)/alpha(ir,it,1)
            IF(tracer_proc(i)/=myproc) THEN
                  WRITE(8888,*) "-----------------------"
                  WRITE(8888,*) "DEBUG FOR INTEGRATION"
                  WRITE(8888,*) "TIME=",time
                  WRITE(8888,*) "Local ID=",i
                  WRITE(8888,*) "V_1=",v_1_t
                  WRITE(8888,*) "V_2=",v_2_t
                  WRITE(8888,*) "Global ID=",tracer_id(i)
                  WRITE(8888,*) "TRACER_PROC =",tracer_proc(i)
                  WRITE(8888,*) "r=",pos_tracer_r(i),tracer_ci_r(i)
                  WRITE(8888,*) "th=",pos_tracer_theta(i),
     &            tracer_ci_theta(i)
                  WRITE(8888,*) "INTERP VALUE FOR THE TRACER V1:"
                  WRITE(8888,*) v_1(ir-1,it-1,1),v_1(ir-1,it,1)
     &             ,v_1(ir,it-1,1),v_1(ir,it,1) 
                  WRITE(8888,*) "INTERP VALUE FOR THE TRACER V2:"
                  WRITE(8888,*) v_2(ir-1,it-1,1),v_2(ir-1,it,1)
     &             ,v_2(ir,it-1,1),v_2(ir,it,1) 
            END IF
            IF(rkint_order.eq.1)    THEN
                  pos_tracer_r_old(i) = pos_tracer_r(i)
                  pos_tracer_r(i) = pos_tracer_r_old(i) + 
     &            v_1_t*delta_t
                  pos_tracer_theta_old(i) = pos_tracer_theta(i)
                  pos_tracer_theta(i) = pos_tracer_theta_old(i) +
     &            v_2_t*delta_t
            ELSE IF(rkint_order.eq.2)     THEN
                  pos_tracer_r(i) = (beta_rk * pos_tracer_r_old(i)
     &            + pos_tracer_r(i)+ v_1_t * delta_t)/alpha_rk
                  pos_tracer_theta(i) = (beta_rk * pos_tracer
     &            _theta_old(i) + pos_tracer_theta(i) + v_2_t
     &            * delta_t)/alpha_rk
            ELSE IF(rkint_order.eq.3)     THEN
                  pos_tracer_r(i) = 
     &            (pos_tracer_r_old(i)+2*(pos_tracer_r(i)
     &            +v_1_t * delta_t)) / 3.0_rk
                  pos_tracer_theta(i) = (pos_tracer_theta_old(i)
     &            + 2*(pos_tracer_theta(i) + v_2_t * delta_t))/
     &            3.0_rk
            END IF
      END DO LOOP_INT
      WRITE(9999,*) "INTEGRATION FINISHED"
      END SUBROUTINE

c     ==================================================================

      SUBROUTINE parameter_interp()

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE configure
      USE gfloat_hy,ONLY: time
      IMPLICIT NONE
c     Currently we're using linear interpolation for every parameter
c     But in further study we may use more accurate interpolation methods like Cubic Spline Interpolation
      
      INTEGER (KIND=ik) :: tci,i,j,tci_t
      REAL (KIND=rk) :: dr_t,dth_t,sr
      LOOP_INTERP:DO i=1, tracer_num
            tci = tracer_ci_r(i)
            tci_t = tracer_ci_theta(i)
            IF (tci.eq.-1) THEN
                  CYCLE LOOP_INTERP
            END IF
            IF ((tci_t.lt.n_s).or.(tci_t.gt.n_e)) THEN
                      CYCLE LOOP_INTERP
            END IF

            dr_t = (pos_tracer_r(i)-r_if(tci-1))/(delta_r(tci))
            dth_t = (pos_tracer_theta(i)-theta_if(tci_t-1))
 
            v_1_tracer(i) = v_1(tci-1,tci_t-1,1) +
     &      (v_1(tci,tci_t-1,1)-v_1(tci-1,tci_t-1,1)) * dr_t +
     &      (v_1(tci-1,tci_t,1)-v_1(tci-1,tci_t-1,1)) * dth_t
            
            v_2_tracer(i) = v_2(tci-1,tci_t-1,1) +
     &      (v_2(tci,tci_t-1,1)-v_2(tci-1,tci_t-1,1)) * dr_t +
     &      (v_2(tci-1,tci_t,1)-v_2(tci-1,tci_t-1,1)) * dth_t
            
            v_1_tracer(i) = g_up_11(tci,tci_t,1)*v_1_tracer(i)
     &      -beta_up_1(tci,tci_t,1)/alpha(tci,tci_t,1)

            v_2_tracer(i) = g_up_22(tci,tci_t,1)*v_2_tracer(i)
     &      -beta_up_2(tci,tci_t,1)/alpha(tci,tci_t,1)
            
            temp_tracer(i) = t(tci-1,tci_t-1,1) + 
     &      (t(tci,tci_t-1,1)-t(tci-1,tci_t-1,1)) * dr_t +
     &      (t(tci-1,tci_t,1)-t(tci-1,tci_t-1,1)) * dth_t


            
            LOOP_XNNU:DO j=1,config%qn
              xnnu_tracer(i,j) = xnnu(tci-1,tci_t-1,1,j)+
     &        (xnnu(tci,tci_t-1,1,j)-xnnu(tci-1,tci_t-1,1,j))*dr_t +
     &        (xnnu(tci-1,tci_t,1,j)-xnnu(tci-1,tci_t-1,1,j))*dth_t
            END DO LOOP_XNNU
            den_tracer(i) = rho(tci-1,tci_t-1,1)+
     &      (rho(tci,tci_t-1,1)-rho(tci-1,tci_t-1,1))*dr_t+
     &      (rho(tci-1,tci_t,1)-rho(tci-1,tci_t-1,1))*dth_t

            pressure_tracer(i) = p(tci-1,tci_t-1,1)
     &      +(p(tci,tci_t-1,1)-p(tci-1,tci_t-1,1))*dr_t
     &      +(p(tci-1,tci_t,1)-p(tci-1,tci_t-1,1))*dth_t

            entropy_tracer(i) = entropy(tci-1,tci_t-1,1)+
     &      (entropy(tci,tci_t-1,1)-entropy(tci-1,tci_t-1,1))*dr_t +
     &      (entropy(tci-1,tci_t,1)-entropy(tci-1,tci_t-1,1))*dth_t
      END DO LOOP_INTERP

      END SUBROUTINE parameter_interp
c     ==================================================================

      SUBROUTINE parameter_interp2()

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE configure
      USE gfloat_hy,ONLY: time
      IMPLICIT NONE
c     Currently we're using linear interpolation for every parameter
c     But in further study we may use more accurate interpolation methods like Cubic Spline Interpolation
      
      INTEGER (KIND=ik) :: ir,it,i,j
      REAL (KIND=rk) :: dr_t,dth_t,sr
      LOGICAL :: inghz
c      WRITE(9999,*) "TEST12"
      LOOP_INTERP:DO i=1, tracer_num
            inghz = .false.
            ir = tracer_ci_r(i)
            it = tracer_ci_theta(i)
c            WRITE(9999,*) "TEST13",i,ir,it
            IF ((ir.lt.0).or.(ir.gt.m)) THEN
                  CYCLE LOOP_INTERP
            END IF
            IF ((it.lt.n_s-1).or.(it.gt.n_e+1)) THEN
                  CYCLE LOOP_INTERP
            END IF
            IF(it.lt.n_s) THEN
               inghz = .true.
               it = n_s
            ELSE IF(it.gt.n_e) THEN
               inghz = .true.
               it = n_e
            END IF
c            WRITE(9999,*) "TEST14"
c            WRITE(8888,*) "----debug for v_1-----"
            v_1_tracer(i) = bilineary_interpolation
     &      (v_1(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i
     &      ,.false.,"v_1")
c            WRITE(8888,*) "----debug for v_2-----"
            v_2_tracer(i) = bilineary_interpolation
     &      (v_2(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i,
     &      .false.,"v_2")

            
            v_1_tracer(i) = g_up_11(ir,it,1)*v_1_tracer(i)
     &      -beta_up_1(ir,it,1)/alpha(ir,it,1)

            v_2_tracer(i) = g_up_22(ir,it,1)*v_2_tracer(i)
     &      -beta_up_2(ir,it,1)/alpha(ir,it,1)
            if(inghz) then
               v_2_tracer(i) = -v_2_tracer(i)
            end if
c            WRITE(8888,*) "----debug for temp_tracer"
            temp_tracer(i) = bilineary_interpolation
     &      (t(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i
     &      ,.true.,"temperature")
            

c            WRITE(8888,*) "----debug for xnnu"
            LOOP_XNNU:DO j=1,config%qn
            xnnu_tracer(i,j) = bilineary_interpolation
     &      (xnnu(-1:m+1,n_s-1:n_e+1,1,j),-1,m+1,n_s-1,n_e+1,1,1,i
     &      ,.false.,"xnnu")
            END DO LOOP_XNNU
c            WRITE(8888,*) "----debug for density----"
            den_tracer(i) = bilineary_interpolation
     &      (rho(-1:m+1,n_s-1:n_e+1,1),-1,m+1,n_s-1,n_e+1,1,1,i
     &      ,.false.,"density")
c            WRITE(8888,*) "----debug for pressure----"
            pressure_tracer(i) = 
     &      bilineary_interpolation(p(-1:m+1,n_s-1:n_e+1,1),
     &      -1,m+1,n_s-1,n_e+1,1,1,i,.false.,"pressure")
c            WRITE(8888,*) "----debug for entropy----"
            entropy_tracer(i) = 
     &      bilineary_interpolation(entropy(-1:m+1,n_s-1:n_e+1,1),-1,
     &      m+1,n_s-1,n_e+1,1,1,i,.false.,"entropy")
c           WRITE(9999,*) "TEST15"
      END DO LOOP_INTERP

      END SUBROUTINE parameter_interp2
c     ==================================================================

      SUBROUTINE tracer_print()

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE gfloat_hy, ONLY: time
      use intgrs_hy, only : nstep
      IMPLICIT NONE
      
      integer (KIND=ik) :: i,tci_r,tci_t
      IF(time-tracer_printtime.le.0.0001_rk) THEN
            return
      END IF
      tracer_printtime = time
      LOOP_PRINT: DO i=1,tracer_num
            tci_r = tracer_ci_r(i)
            tci_t = tracer_ci_theta(i)
            WRITE(7777,*) "-----------------------"
            WRITE(7777,*) 'TRACER NO = ', i,'TIME = ', time," ",nstep
            WRITE(7777,*) 'GLOBAL ID=',tracer_id(i)
            WRITE(7777,*) 'IPROC = ', tracer_proc(i)
            WRITE(7777,*) 'position r = ',pos_tracer_r(i),tci_r
            WRITE(7777,*) 'position theta = ', pos_tracer_theta(i),tci_t
c            WRITE(7777,*) 'v = ',v_1_tracer(i)
c            WRITE(7777,*) 'v_theta = ', v_2_tracer(i)
c            WRITE(7777,*) 't = ',temp_tracer(i)
c            WRITE(7777,*) 'p = ',pressure_tracer(i)
c            WRITE(7777,*) 's = ',entropy_tracer(i)
c            WRITE(7777,*) 'den = ',den_tracer(i) 
      END DO LOOP_PRINT
c      WRITE(9999,*) "TEST17"
      END SUBROUTINE tracer_print
c     ==================================================================
      FUNCTION bilineary_interpolation(arr,l1,u1,l2,u2,l3,u3,ti,negc
     & ,dbg) result(val)

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE gfloat_hy, ONLY: time
      INTEGER, INTENT(in) :: ti,l1,u1,l2,u2,l3,u3
      REAL  :: arr(l1:u1,l2:u2,l3:u3)
      REAL :: val
      REAL :: xir0,xir1,xit0,xit1,pos_r,pos_t
      INTEGER :: ir,it
      LOGICAL, INTENT(in) :: negc
      CHARACTER(len=*), INTENT(in) :: dbg
      

      pos_r = pos_tracer_r(ti)
      pos_t = pos_tracer_theta(ti)
      ir = tracer_ci_r(ti)
      it = tracer_ci_theta(ti)
      xir0 = (pos_r-r_if(ir-1))/delta_r(ir)
      xir1 = 1 - xir0
      xit0 = (pos_t-theta_if(it-1))/delta_theta
      xit1 = 1 - xit0
      if(xit0.lt.0) then
            xit0 = xit0-floor(xit0)
            xit1 = 1-xit0
      else  if(xit1.lt.0) then
            xit1 = xit1-floor(xit1)
            xit0=1-xit1
      end if
      val =  xit1*(xir1*arr(ir-1,it-1,1)+xir0*arr(ir,it-1,1))+xit0*
     &     (xir1*arr(ir-1,it,1)+xir0*arr(ir,it,1))
      IF(negc.and.neg_check(ti)) THEN
            IF(val.lt.-1e-20) THEN
      neg_check(ti) = .false.
      WRITE(8888,*) dbg
      WRITE(8888,*) "basic parameters:",ti
      WRITE(8888,*) "l1:",l1,"l2:",l2,"l3:",l3
      WRITE(8888,*) "u1:",u1,"u2:",u2,"u3:",u3
      WRITE(8888,*) "pos_r:",pos_r,"pos_t:",pos_t
      WRITE(8888,*) "delta_r",delta_r(ir),"delta_t",delta_theta
      WRITE(8888,*) "ir:",ir,"it:",it
      WRITE(8888,*) "t_if",theta_if(it-1)
      WRITE(8888,*) "t_ir",r_if(ir-1)
      WRITE(8888,*) "xir0:",xir0,"xir1:",xir1
      WRITE(8888,*) "xit0:",xit0,"xit1:",xit1
      WRITE(8888,*) "biinterp values:"
      WRITE(8888,*) arr(ir-1,it-1,1)
      WRITE(8888,*) arr(ir-1,it,1)
      WRITE(8888,*) arr(ir,it-1,1)
      WRITE(8888,*) arr(ir,it,1)
      WRITE(8888,*) val                  
            END IF
      END IF
      IF(nan_check(ti)) THEN
            IF(isnan(val)) THEN  
      nan_check(ti) = .false.
      WRITE(8888,*) dbg
      WRITE(8888,*) "basic parameters:",ti
      WRITE(8888,*) "l1:",l1,"l2:",l2,"l3:",l3
      WRITE(8888,*) "u1:",u1,"u2:",u2,"u3:",u3
      WRITE(8888,*) "pos_r:",pos_r,"pos_t:",pos_t
      WRITE(8888,*) "delta_r",delta_r(ir),"delta_t",delta_theta
      WRITE(8888,*) "ir:",ir,"it:",it
      WRITE(8888,*) "t_if",theta_if(it-1)
      WRITE(8888,*) "t_ir",r_if(ir-1)
      WRITE(8888,*) "xir0:",xir0,"xir1:",xir1
      WRITE(8888,*) "xit0:",xit0,"xit1:",xit1
      WRITE(8888,*) "biinterp values:"
      WRITE(8888,*) arr(ir-1,it-1,1)
      WRITE(8888,*) arr(ir-1,it,1)
      WRITE(8888,*) arr(ir,it-1,1)
      WRITE(8888,*) arr(ir,it,1)
      WRITE(8888,*) val                              
            END IF
      END IF
      END FUNCTION bilineary_interpolation 
c     ==================================================================
      SUBROUTINE MPI_TRACER()

c     ==================================================================
      USE precision
      USE grid_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE tracer_cfc
      USE gfloat_hy, ONLY: time
      USE mo_mpi
      use phycon, only : pc_pi
      use intgrs_hy, only : nstep
      INTEGER :: istat
      INTEGER(kind=ik) :: nsend_l,nsend_r,nrecv_l,nrecv_r,tracer_num_new
      INTEGER(kind=ik) :: nsl,nsr,nns,ndebug
      INTEGER(kind=ik) :: i,j,k
C      LOGICAL :: issent,numchanged
c     BM: More efficient to use an array for issent
      LOGICAL :: numchanged
      LOGICAL :: issent (1:tracer_num)
      REAL(kind=rk), allocatable :: sbuffer_tl(:,:)
      REAL(kind=rk), allocatable :: rbuffer_tl(:,:)
      REAL(kind=rk), allocatable :: sbuffer_tr(:,:)
      REAL(kind=rk), allocatable :: rbuffer_tr(:,:)
      
      REAL(kind=rk), allocatable :: buffer_tnew(:,:)

      INTEGER (kind=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      
      integer (kind=ik), parameter :: tag_ysndl=9700_ik
      integer (kind=ik), parameter :: tag_ysndr=9701_ik
      integer (kind=ik), parameter :: tag_zsndl=9702_ik
      integer (kind=ik), parameter :: tag_zsndr=9703_ik
      nsend_l = 0
      nsend_r = 0
      numchanged = .false.
c     BM: Initially set issent to false.
      issent (:) = .false.
      WRITE(9999,*) "START MPI TRACER TIME = ",time," ",nstep
      do i=1,tracer_num
            call tci_search(i)
            call tci_t_search(i)
            if(pos_tracer_theta(i).lt.theta_if(n_s-1)
     &      .and.pos_tracer_theta(i).gt.theta_if(0)) then
                  nsend_l = nsend_l + 1
            end if
            if(pos_tracer_theta(i).gt.theta_if(n_e)
     &      .and.pos_tracer_theta(i).lt.theta_if(n)) then
                  nsend_r = nsend_r + 1
            end if
      end do
      WRITE(9999,*) "NSEND_L:",nsend_l,"NSEND_R:",nsend_r
      allocate(sbuffer_tl(1:nsend_l,1:5),sbuffer_tr(1:nsend_r,1:5),
     &  stat=istat)
      WRITE(9999,*) "COUNTING SEND FINISHED"
      nsl=0
      nsr=0
      do i = 1,tracer_num
C     if(pos_tracer_theta(i).ge.theta_if(0)) then
c     BM: Change .ge. to .gt. here:
        if(pos_tracer_theta(i).gt.theta_if(0)) then
          if(pos_tracer_theta(i).lt.theta_if(n_s-1)) then
            nsl = nsl+1
            sbuffer_tl(nsl,1) = tracer_id(i)
            sbuffer_tl(nsl,2) = pos_tracer_r(i)
            sbuffer_tl(nsl,3) = pos_tracer_theta(i)
            sbuffer_tl(nsl,4) = tracer_ci_r(i)
            sbuffer_tl(nsl,5) = tracer_proc(i)
            WRITE(9999,*) "-----------------------"
            WRITE(9999,*) "TRACER TO BE SENT TO LEFT:"
            WRITE(9999,*) "TRACER ID:",tracer_id(i)
            WRITE(9999,*) "POS_R:", pos_tracer_r(i)
            WRITE(9999,*) "POS_THETA:",pos_tracer_theta(i)
            WRITE(9999,*) "TRACER_CI_R",tracer_ci_r(i)
c     BM
            issent (i) = .TRUE.
          end if
        else
          pos_tracer_theta(i) = - pos_tracer_theta(i)
        end if
c        if(pos_tracer_theta(i).ge.theta_if(0)) then
c     BM: I believe the previous line should instead be
        if(pos_tracer_theta(i).lt.theta_if(n)) then
          if(pos_tracer_theta(i).gt.theta_if(n_e)) then
            nsr = nsr + 1
            sbuffer_tr(nsr,1) = tracer_id(i)
            sbuffer_tr(nsr,2) = pos_tracer_r(i)
            sbuffer_tr(nsr,3) = pos_tracer_theta(i)
            sbuffer_tr(nsr,4) = tracer_ci_r(i)
            sbuffer_tr(nsr,5) = tracer_proc(i)
            WRITE(9999,*) "-----------------------"
            WRITE(9999,*) "TRACER TO BE SENT TO RIGHT:"
            WRITE(9999,*) "TRACER ID:",tracer_id(i)
            WRITE(9999,*) "POS_R:", pos_tracer_r(i)
            WRITE(9999,*) "POS_THETA:",pos_tracer_theta(i)
            WRITE(9999,*) "TRACER_CI_R",tracer_ci_r(i)
c     BM
            issent (i) = .TRUE.
          end if
        else
            pos_tracer_theta(i) = 2.0_rk * pc_pi - pos_tracer_theta(i)
        end if
      end do
      if(nsl/=nsend_l.or.nsr/=nsend_r) then
            WRITE(9999,*) "ERROR:SNED NUM IS NOT AS EXPECTED"
      end if
      WRITE(9999,*) "SBUFFER FINISHED"

      
      call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
      call MPI_Sendrecv (nsend_l,1,
     &  MPI_INTEGER,dest,tag_ysndl,nrecv_r,
     &  1,MPI_INTEGER,
     &  src,tag_ysndl,cart_comm,mpistat,ierr)
        if(src.eq.MPI_PROC_NULL) then
          nrecv_r = 0
        end if

      call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
      call MPI_Sendrecv (nsend_r,1,
     &  MPI_INTEGER,dest,tag_ysndr,nrecv_l,
     &  1,MPI_INTEGER,
     &  src,tag_ysndr,cart_comm,mpistat,ierr)

        if(src.eq.MPI_PROC_NULL) then
          nrecv_l = 0
        end if
      WRITE(9999,*) "NUMBER COMMUNICATION FINISHED"
        tracer_num_new = tracer_num + (nrecv_l+nrecv_r) 
     &  - (nsend_l+nsend_r)

      IF(nsend_l/=0.or.nsend_r/=0.or.nrecv_l/=0.or.nrecv_r/=0) THEN
            numchanged = .true.
            WRITE(9999,*) "NSEND_L:",nsend_l,"NSEND_R:",nsend_r
            WRITE(9999,*) "NRECV_L:",nrecv_l,"NRECV_R:",nrecv_r
            WRITE(9999,*) "TRACER_NUM_NEW:",tracer_num_new
            WRITE(9999,*) "TRACER_NUM",tracer_num
      END IF
c      WRITE(9999,*) "TEST4"
      allocate(buffer_tnew(1:tracer_num_new,1:5),
     &  stat=istat)
c     WRITE(9999,*) "TEST5"
      nns = 0
      ndebug = 0
      LOOP_Main:do i = 1,tracer_num
c     BM: You can code this more efficiently
c$$$        issent = .false.
c$$$        LOOP_FINDL:do j=1,nsend_l 
c$$$            if(tracer_id(i).eq.sbuffer_tl(j,1)) then
c$$$                issent = .true.
c$$$                EXIT LOOP_FINDL
c$$$            end if
c$$$        end do LOOP_FINDL
c$$$        LOOP_FINDR: do j=1,nsend_r
c$$$            if(tracer_id(i).eq.sbuffer_tr(j,1)) then
c$$$                issent = .true.
c$$$                EXIT LOOP_FINDR
c$$$            end if
c$$$        end do LOOP_FINDR
c$$$        if(.not.issent) then
        if(.not.issent(i)) then
            nns = nns+1
            buffer_tnew(nns,1) = tracer_id(i)
            buffer_tnew(nns,2) = pos_tracer_r(i)
            buffer_tnew(nns,3) = pos_tracer_theta(i)
            buffer_tnew(nns,4) = tracer_ci_r(i)
            buffer_tnew(nns,5) = tracer_proc(i)
        end if
      end do LOOP_Main

      if(nns /= (tracer_num-(nsend_l+nsend_r))) then
        WRITE(9999,*) "ERROR: NUMBER OF TRACERS NOT SENT IS NOT AS 
     &  EXPECTED",nns,tracer_num-(nsend_l+nsend_r),time
      end if

      allocate(rbuffer_tl(1:nrecv_l,1:5),
     &  rbuffer_tr(1:nrecv_r,1:5),stat=istat)

      call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
      call MPI_Sendrecv (sbuffer_tl,nsend_l*5,
     &  MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbuffer_tr,
     &  nrecv_r*5,MPI_DOUBLE_PRECISION,
     &  src,tag_ysndl,cart_comm,mpistat,ierr)

      call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
      call MPI_Sendrecv (sbuffer_tr,nsend_r*5,
     &  MPI_DOUBLE_PRECISION,dest,tag_ysndr,rbuffer_tl,
     &  nrecv_l*5,MPI_DOUBLE_PRECISION,
     &  src,tag_ysndr,cart_comm,mpistat,ierr)

      WRITE(9999,*) "BUFFER COMMUNICATION FINISHED"

      if((nns + nrecv_l +nrecv_r)/=tracer_num_new) then
        WRITE(9999,*) "ERROR: NUMBER OF TRACERS RECEIVED IS NOT AS 
     &  EXPECTED",nns,nrecv_l,nrecv_r,tracer_num_new,time
      end if
      if(nrecv_l/=0) then
            WRITE(9999,*) "WRITE RECEIVED FROM LEFT"
            do i=1,nrecv_l
                  WRITE(9999,*) rbuffer_tl(i,1)
                  WRITE(9999,*) rbuffer_tl(i,2)
                  WRITE(9999,*) rbuffer_tl(i,3)
                  WRITE(9999,*) rbuffer_tl(i,4)
                  WRITE(9999,*) rbuffer_tl(i,5)
            end do
      end if
      if(nrecv_r/=0) then
            WRITE(9999,*) "WRITE RECEIVED FROM RIGHT"
            do i=1,nrecv_r
                  WRITE(9999,*) rbuffer_tr(i,1)
                  WRITE(9999,*) rbuffer_tr(i,2)
                  WRITE(9999,*) rbuffer_tr(i,3)
                  WRITE(9999,*) rbuffer_tr(i,4)
                  WRITE(9999,*) rbuffer_tr(i,5)
            end do
      end if
      buffer_tnew(nns+1:nns+nrecv_l,:) = rbuffer_tl(:,:)
      nns = nns+  nrecv_l
      buffer_tnew(nns+1:tracer_num_new,:) = rbuffer_tr(:,:)

      deallocate(tracer_id,pos_tracer_r,pos_tracer_theta,tracer_ci_r
     & ,tracer_proc)

      allocate(
     &  tracer_id(tracer_num_new),
     &  pos_tracer_r(tracer_num_new),
     &  pos_tracer_theta(tracer_num_new),
     &  tracer_ci_r(tracer_num_new),
     &  tracer_proc(tracer_num_new))

      tracer_id = buffer_tnew(:,1)
      pos_tracer_r = buffer_tnew(:,2)
      pos_tracer_theta = buffer_tnew(:,3)
      tracer_ci_r = buffer_tnew(:,4)
      tracer_proc = buffer_tnew(:,5)
      tracer_num = tracer_num_new
      WRITE(9999,*) "New Tracers NUM:",tracer_num_new
      deallocate(sbuffer_tl,sbuffer_tr,rbuffer_tl,
     &  rbuffer_tr,buffer_tnew)
      do i=1,tracer_num_new
            WRITE(9999,*) "-----------------------"
            WRITE(9999,*) "local id:",i
            WRITE(9999,*) "global_id",tracer_id(i)
            WRITE(9999,*) "pos_r:",pos_tracer_r(i)
            WRITE(9999,*) "pos_theta:",pos_tracer_theta(i)
            WRITE(9999,*) "radial ci index:",tracer_ci_r(i)
            WRITE(9999,*) "FROM PROC:",tracer_proc(i)
      end do
      END SUBROUTINE MPI_TRACER
c  ********************************************************************
      END MODULE tracer
