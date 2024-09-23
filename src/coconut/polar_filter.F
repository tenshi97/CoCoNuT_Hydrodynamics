c     ==================================================================

      MODULE filtering

c     ==================================================================
      
      CONTAINS


c     ==================================================================

      SUBROUTINE polar_filter (mode)

c     ==================================================================
c     - Filter fluxes and pressure on varphi interfaces.
c     ------------------------------------------------------------------

      USE precision
      
      USE conserved_cfc
      USE fluxes_cfc
      USE grid_cfc
      USE hydro_areas_mod
      USE hydro_primitives_cfc, ONLY: p
      USE perm_aux_cfc, ONLY: delta_theta_inv, sinus_theta,
     &     sinus_theta_if
      USE phycon, ONLY: pc_pi
      USE size_cfc

      USE fft_mpi
      USE mo_mpi

      USE configure

      USE, INTRINSIC :: iso_c_binding

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mode

      REAL (KIND=rk) :: xin (1:m*(6+config%qn), n_s:n_e, o_s:o_e)
      REAL (KIND=rk) :: xout (1:m*(6+config%qn), n_s:n_e, o_s:o_e)

      REAL (KIND=rk) :: sbufz (1:m, n_s:n_e, 6+config%qn)
      REAL (KIND=rk) :: rbufz (1:m, n_s:n_e, 6+config%qn)

      INTEGER (KIND=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      INTEGER (KIND=ik), parameter :: tag_ysndl=1500_ik
      INTEGER (KIND=ik), parameter :: tag_ysndr=1501_ik
      INTEGER (KIND=ik), parameter :: tag_zsndl=1502_ik
      INTEGER (KIND=ik), parameter :: tag_zsndr=1503_ik


      REAL (KIND=rk) :: wn, rspher, scr, sth
      REAL (KIND=rk) :: dsmooth (1:m)
      INTEGER (KIND=ik) :: index (1:o)

      INTEGER (KIND=ik) :: i, ii, j, k, jk, l, ixf, ioy
      INTEGER (KIND=ik) :: msmooth (n_s:n_e)

#if DIMN==3
      
      ixf    = areas%ix_are(1, 2)
      ioy    = areas%ix_are(1, 6)
      if (ioy.eq.n) then
         rspher = r (ixf + 1)
      else
         rspher = 1e6
      end if

      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         IF (mode .EQ. 0) THEN
            DO l = 1, 5
               xin (1+(l-1)*m:l*m,j,k)=f_num_if(1:m,j,k,l)
            END DO
            DO l = 1, config%qn
               xin (1+(l+4)*m:(l+5)*m,j,k)=f_num_xnu_if(1:m,j,k,l)
            END DO
            xin (1+(config%qn+5)*m:(config%qn+6)*m,j,k)=
     &           pav_varphi(1:m,j,k)
         ELSE
            xin (1+0*m:1*m,j,k)=d_cap_hat(1:m,j,k)
            xin (1+1*m:2*m,j,k)=s_1_hat  (1:m,j,k)
            xin (1+2*m:3*m,j,k)=s_2_hat  (1:m,j,k)
            xin (1+3*m:4*m,j,k)=s_3_hat  (1:m,j,k)
            xin (1+4*m:5*m,j,k)=tau_hat  (1:m,j,k)
            DO l = 1, config%qn
               xin (1+(l+4)*m:(l+5)*m,j,k)=d_cap_xnu_hat(1:m,j,k,l)
            END DO
            xin (1+(config%qn+5)*m:(config%qn+6)*m,j,k)=
     &           p(1:m,j,k)
         END IF

      END DO

      CALL FFT (xin (1:m*(6+config%qn), n_s:n_e, o_s:o_e),
     &     xout (1:m*(6+config%qn), n_s:n_e, o_s:o_e), m*(6+config%qn), 
     &     1, index)
         
      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         wn = REAL (MIN (index (k), o - index (k)), KIND=rk)
         wn = MAX (wn, 1e-35_rk)
         sth = MIN (sinus_theta_if (j-1), sinus_theta_if (j)) 
         sth = MAX (sth, 0.5_rk*sinus_theta (1))

         IF (mode .LE. 1) THEN
            msmooth (j) = m
            DO i = 1, m
               IF (sth * delta_varphi * 
     &              delta_theta_inv *
     &              r (i) / rspher * 2.0_rk *
     &              MAX(0.7_rk / config%cfl, 1.0_rk) .LT. 1.0_rk)
     &              msmooth (j) = i
            END DO
c            PRINT *,'msmooth',j,msmooth(j)

            DO ii = 0, (5+config%qn)*m, m
               DO i = 1, m
                  xin (i+ii, j, k) =
     &                 xout (i+ii, j, k) *
     &                 MIN (sth * delta_varphi * 
     &                 delta_theta_inv *
     &                 r (i) / rspher /
     &                 SIN (0.5_rk * wn * delta_varphi) * 2.0_rk *
     &                 MAX(0.7_rk / config%cfl, 1.0_rk),
     &                 1.0_rk) ** 2
               END DO
               IF ((j.EQ.1 .OR. j.EQ.n) .AND. wn.GE.0.9_rk)
     &              xin (ii+1:ii+msmooth(j), j, k) = 0.0_rk
            END DO
 
         ELSE
#if 1
            dsmooth(:) = sth*delta_varphi*delta_theta_inv*0.5_rk
     &           *r(1:m)/rspher
            DO i = 1, m
               IF (dsmooth(i) .LE. 0.5_rk) THEN
c                  dsmooth(i) = 0.25_rk * delta_varphi / dsmooth
                  dsmooth(i) = 0.2_rk * delta_varphi / dsmooth(i)
               ELSE
                  dsmooth(i) = 0.0_rk
               ENDIF
            END DO
            IF (wn .GT. 0.9_rk .AND. (j .EQ. 1 .OR. j .EQ.n)) THEN
               DO ii = 0, (5+config%qn)*m, m
                  DO i = 1, m !smooth (j)
                     xin (i + ii, j, k) = 0.0_rk
                  END DO
               END DO
            ELSE
               DO ii = 0, (5+config%qn)*m, m
                  DO i = 1, m
                     xin (i + ii, j, k) =
     &                    xout (i + ii, j, k) *
     &                    ((SIN(0.5_rk * dsmooth(i) * wn) + 1e-40_rk) / 
     &                    (0.5_rk * dsmooth(i) * wn + 1e-40_rk)) ** 3
                  END DO
               END DO
            ENDIF
#else
            msmooth (j) = m
            DO i = 1, m
               IF (sth * delta_varphi * 
     &              delta_theta_inv *
     &              r (i) / rspher * 2.0_rk *
     &              MAX(0.7_rk / config%cfl, 1.0_rk) .LT. 1.0_rk)
     &              msmooth (j) = i
            END DO
            PRINT *,'msmooth',j,msmooth(j)

            DO ii = 0, (5+config%qn)*m, m
               DO i = 1, m
                  xin (i+ii, j, k) =
     &                 xout (i+ii, j, k) *
     &                 MIN (sth * delta_varphi * 
     &                 delta_theta_inv *
     &                 r (i) / rspher /
     &                 SIN (0.5_rk * wn * delta_varphi) * 2.0_rk *
     &                 MAX(0.7_rk / config%cfl, 1.0_rk),
     &                 1.0_rk) ** 2
               END DO
               IF ((j.EQ.1 .OR. j.EQ.n) .AND. wn.GE.0.9_rk)
     &              xin (ii+1:ii+msmooth(j), j, k) = 0.0_rk
            END DO
#endif
         END IF

      END DO

      CALL FFT (xin (1:m*(6+config%qn), n_s:n_e, o_s:o_e),
     &     xout (1:m*(6+config%qn), n_s:n_e, o_s:o_e), m*(6+config%qn),
     &     -1, index)
      
      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         IF (mode .EQ. 0) THEN
            DO l = 1, 5
               f_num_if(1:m,j,k,l) = xout (1+(l-1)*m:l*m,j,k)
            END DO
            DO l = 1, config%qn
               f_num_xnu_if(1:m,j,k,l) = xout (1+(l+4)*m:(l+5)*m,j,k)
            END DO
            pav_varphi(1:m,j,k) = 
     &           xout (1+(config%qn+5)*m:(config%qn+6)*m,j,k)
         ELSE
            d_cap_hat(1:m,j,k)=xout (1+0*m:1*m,j,k)
            s_1_hat  (1:m,j,k)=xout (1+1*m:2*m,j,k)
            s_2_hat  (1:m,j,k)=xout (1+2*m:3*m,j,k)
            s_3_hat  (1:m,j,k)=xout (1+3*m:4*m,j,k)
            tau_hat  (1:m,j,k)=xout (1+4*m:5*m,j,k)
            DO l = 1, config%qn
               d_cap_xnu_hat(1:m,j,k,l)=xout (1+(l+4)*m:(l+5)*m,j,k)
            END DO
            IF (mode .EQ. 2) THEN
            p(1:m,j,k) =
     &       xout (1+(config%qn+5)*m:(config%qn+6)*m,j,k)
         END IF
         END IF
      END DO
         


      IF (mode .EQ. 0) THEN
c     Communicate required fluxes to right neighbour
         DO l = 1, 5
            sbufz(1:m,n_s:n_e,l)=
     &        f_num_if (1:m,n_s:n_e,o_e,l)
         END DO
         DO l = 1, config%qn
            sbufz(1:m,n_s:n_e,l+5)=
     &           f_num_if (1:m,n_s:n_e,o_e,l)
         END DO
         sbufz(1:m,n_s:n_e,config%qn+6) = pav_varphi (1:m,n_s:n_e,o_s-1)

c     determine source and destination ranks and send data
         CALL MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
            
         IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
            CALL MPI_Send (sbufz,m*n_loc*(6_ik+config%qn),
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
            CALL MPI_Recv (rbufz,m*n_loc*(6_ik+config%qn),
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
         ELSE
            CALL MPI_Recv (rbufz,m*n_loc*(6_ik+config%qn),
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
            CALL MPI_Send (sbufz,m*n_loc*(6_ik+config%qn),
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
         END IF
         
c      PRINT *,'fft',k,f_num_if(285,n_s,o_s-1:o_e,1),
c     &     xout (285,n_s,o_s:o_e)
         DO l = 1, 5
            f_num_if (1:m,n_s:n_e,o_s-1,l) =
     &           rbufz (1:m,n_s:n_e,l)
         END DO
         DO l = 1, config%qn
            f_num_if (1:m,n_s:n_e,o_s-1,l) =
     &           rbufz (1:m,n_s:n_e,l+5)
         END DO
         pav_varphi (1:m,n_s:n_e,o_s-1) = rbufz(1:m,n_s:n_e,config%qn+6)
      END IF
      
C      STOP
#endif /* DIMN==3 */
      
      RETURN

      END SUBROUTINE polar_filter

c     ==================================================================



      END MODULE filtering

c     ==================================================================



