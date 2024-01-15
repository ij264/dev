      program driver
      use constants
      use ref_flds, only : ref_flds_init
      use sph_tools, only : init_sph_tls, theta, phi
      use util, only : r, rho0, rvsc, initialize, viscinit
      use commpi
      implicit none

      common /ResCorr/ GeoCorr,GravCorr,DyntCorr,CMBCorr, surmax, geomax
      common /path/     output, input, name_prof

      interface 
         subroutine inviscpll
         end subroutine inviscpll
      end interface

      integer N, sufi, suff, step
      real(dp) time_1, time_2
      real(dp) GeoCorr, GravCorr, DyntCorr, CMBCorr, surmax, geomax
      real(dp) start, finish
      character*120 output, input, name_prof


      N = IARGC()

      call cpu_time(time_1)
      if (N .ne. 3 .and. 
     &   (.not. ((N .eq.2) .and. (inv_flg==1)))) then
         write(*,'(A)') 'Number of given arguments should be 3!'
         write(*,'(I1.1, A)') N, ' is given. Try the following order:'
         write(*,'(A)') './GEOID [INPUT] [OUTPUT] [VISC_PROFILE]'
         stop
      endif

      call getarg(1,input)       ! input is the density file
      call getarg(2,output)      ! ouput is the path to output tomographic_models
      call getarg(3,name_prof)   ! viscosity profile

      call cpu_time(start)
      
      call system('mkdir -p '//trim(output))

      ! initialize the parallel process
      call pinit

      !> reading the density and redius
      call initialize
      !> initialise the reference fields
      call ref_flds_init
      !> initialise the legendre polynomials
      call init_sph_tls 
      !> initialize the necessary parameters and constants, in module consts
      call init_const
      !> we decide if are doing a Monte-Carlo inversion or not 
      if (inv_flg == 1) then
         call monte_carlo
      else if (inv_flg == 2) then
         call viscinit
         call inviscpll 
      else
         call viscinit
         !> compute kernels using propagater matrix tech.
         call prop_matrix
         !> write out the computed correlation
         if (CorrFlag == 1 .and. myproc==root_proc) then
            write(*,'(A)') '-------------------------------------'
            write(*,'(A, I2.1)') ' Processor reporting:', myproc
            write(*,'(A)') '-------------------------------------'
            write(*,'( 2X, A,/, 27X, f15.3)')
     &        'Geoid Correlation:',GeoCorr
            write(*,'( 2X, A,/, 27X, f15.3)')
     &        'Free-Air Gravity Correlation:', GravCorr
            write(*,'( 2X, A,/, 27X, f15.3)')
     &        'Dynamic Topography Correlation:', DyntCorr
            write(*,'( 2X, A,/, 27X, f15.3)')
     &        'CMB     Topography Correlation:', CMBCorr
            write(*,'(2X, A,/, 27X, f15.3)')
     &        'Max Geoid:',geomax
            write(*,'(2X, A,/, 27X, f15.3 )')
     &        'Max SurfTopography:',surmax
            write(*,'(A)') '-------------------------------------'
            call cpu_time(finish)
         endif
      endif

      call cpu_time(time_2)
      if (myproc==root_proc) write(6,'(A,I2.2,A,I2.2,A,I2.2,A)')
     & '  Elapsed Time:                 ',
     & int((time_2-time_1)/3600), ':',
     &   int((time_2 - time_1)/60),''':',mod(int(time_2-time_1),60),'"'

      call pexit
      end program driver

*dk subroutine monte_carlo
      subroutine monte_carlo
      
      use constants
      use util, only : linspace, matrix_out, set_simple_visc
      implicit none
      
      common /path/     output, input, name_prof
      common /ResCorr/ GeoCorr,GravCorr, DyntCorr,CMBCorr, surmax,geomax
      common /ResStat/ flatness, obs_surmax, obs_gemax
      
      character*120 output, input, name_prof
      real(dp)    GeoCorr, GravCorr, DyntCorr,CMBCorr, surmax, geomax
      real(dp)    obs_surmax, obs_gemax
      real(dp)    time_1, time_2, flatness
      real(qp)    all_visc(int_visc), all_thck(int_thck)
      real(qp)    geoid_corr_arr (int_visc, int_thck)
      real(qp)    grav_corr_arr  (int_visc, int_thck)
      real(qp)    dyntop_corr_arr(int_visc, int_thck)
      real(qp)    surmax_arr     (int_visc, int_thck)
      real(qp)    geomax_arr     (int_visc, int_thck)
      real(qp)    flatness_arr   (int_visc, int_thck)
      integer(i4b) ir, ii, ir_lith

      geoid_corr_arr (:,:) = 0.0 
      grav_corr_arr  (:,:) = 0.0
      dyntop_corr_arr(:,:) = 0.0
      surmax_arr     (:,:) = 0.0
      geomax_arr     (:,:) = 0.0

      call linspace(all_visc, 
     &      log10(min_visc_asth), log10(visc_mantle), int_visc)
      all_visc = 10**(all_visc)
      open(001, file=trim(output)//'/ALL_visc',
     &      action='write')
      call matrix_out(all_visc, int_visc, 1, 001)
      close(001)

      call linspace(all_thck, 
     &      min_thck_asth, max_thck_asth, int_thck)
      open(001, file=trim(output)//'/ALL_thck',
     &      action='write')
      call matrix_out(all_thck, 1, int_thck, 001)
      close(001)

      ir_lith = minloc(abs(r(:) - (r(nr+1) - thck_lith)),1) 

      call cpu_time(time_1)
      do ii=1, int_visc
         write(*,'(A,I2.2, A, I2.2)')
     &         'Iteration: ', ii,'/', int_visc
         GeoCorr  = 0.0
         GravCorr = 0.0
         DyntCorr = 0.0
         surmax   = 0.0
         geomax   = 0.0
         flatness = 0.0
         do ir=1,int_thck
            call set_simple_visc(all_visc(ii), all_thck(ir), ir_lith)
            call prop_matrix
            geoid_corr_arr(ii, ir)  = real(GeoCorr    , kind=qp) 
            grav_corr_arr(ii, ir)   = real(GravCorr   , kind=qp) 
            dyntop_corr_arr(ii, ir) = real(DyntCorr   , kind=qp)
            surmax_arr(ii, ir)      = real(surmax     , kind=qp)
            geomax_arr(ii, ir)      = real(geomax     , kind=qp)
            flatness_arr(ii, ir)    = real(flatness   , kind=qp)
         enddo
      enddo
      call cpu_time(time_2)
      write(*,'(A, f6.2, A)') 
     & 'Total elapsed time: ',(time_2 - time_1)/60, ' [min]'

      open(001, file=trim(output)//'/CORR_Geoid_corr.dat',
     &      action='write')
      open(002, file=trim(output)//'/CORR_Grav_corr.dat',
     &      action='write')
      open(003, file=trim(output)//'/CORR_DynT_corr.dat',
     &      action='write')
      open(004, file=trim(output)//'/MAX_surmax.dat',
     &      action='write')
      open(005, file=trim(output)//'/MAX_geomax.dat',
     &      action='write')
      open(006, file=trim(output)//'/FLATNESS.dat',
     &      action='write')

      call matrix_out(geoid_corr_arr,  int_visc, int_thck, 001)
      call matrix_out(grav_corr_arr,   int_visc, int_thck, 002)
      call matrix_out(dyntop_corr_arr, int_visc, int_thck, 003)
      call matrix_out(surmax_arr,      int_visc, int_thck, 004)
      call matrix_out(geomax_arr,      int_visc, int_thck, 005)
      call matrix_out(flatness_arr ,   int_visc, int_thck, 006)

      close(001)
      close(002)
      close(003)
      close(004)
      close(005)
      close(006)

      end subroutine monte_carlo

*dk subroutine prop_matrix
      subroutine prop_matrix

      use commpi
      use constants
      use util, only : r, rho0, rvsc, kernel_out, summary
      implicit none
      common /path/ output, input, name_prof
      integer i, j, k, m, n, ierr, iiir
      integer ipiv(8), kki, kkj
      character*120 output, input, name_prof
      real(dp) GeoKernel(l_min:l_max_calc,nr+1)
      real(dp) VelKernel(l_min:l_max_calc,nr+1)
      real(dp) SurKernel(l_min:l_max_calc,nr+1)
      real(dp) GrvKernel(l_min:l_max_calc,nr+1)
      real(dp) ADMKernel(l_min:l_max_calc,nr+1)
      real(dp) CMBKernel(l_min:l_max_calc,nr+1)
      real(qp) kappa
      real(qp) a(6,6), prop(6,6,1:nr+2)
      real(qp) eigen_i(6), eigen_r(6)
      real(qp) c(8,8), b(8), u(8)
      real(qp) dump, fac
      real(qp) grav
      real(qp) rcond 

      GeoKernel(:,:) = 0.0
      VelKernel(:,:) = 0.0
      CMBKernel(:,:) = 0.0
      SurKernel(:,:) = 0.0
      GrvKernel(:,:) = 0.0
      ADMKernel(:,:) = 0.0
!      if (myproc==root_proc) then
!         print *, rvsc(5)
!      endif
      !> constructing the matrix
      do j= l_min, l_max_calc, 1
      if (mod(j-l_min,nproc) == myproc) then
         !> initialisation of a 6*6 matrix with zero values
         call matunitinit(prop(:,:,nr+2), 6,6)
         do i=nr+1,2,-1
            call aconstruct(a, i, real(j,kind=qp))
            call eigen(a, 6, eigen_r, eigen_i)
            call propmat(a, prop(:,:,i), eigen_r,eigen_i,r(i)/r(i-1))
            prop(:,:,i) = matmul(prop(:,:,i+1), prop(:,:,i))
         enddo
         !> prop(:,:,2) represents prop{6*6}(c,a)
         call cconstruct(c, j, prop(:,:,2))
         call qgeco(c,8,8, ipiv, rcond, u)
         if (rcond .lt. 0.0) then
            write(*,*) rcond
            stop
         endif
         do i=nr+1,1,-1
            !call bconstruct(i, prop(:,:,i+1), j, b)
            call bconstruct(i, prop(:,:,i+1), j, u)
            !call solve_lgs(c, b, 8, u, ipiv)
            call qgesl(c, 8, 8, ipiv, u,0)
            !> geoid kernel
            GeoKernel(j,i)       =
     &        (visc0/(bgrho*r(nr+1)*grav(r(nr+1))))*u(3)
            !> free air gravity anomaly kernel, (in our idealized geometry
            !>    gravity anomaly is the same as free air gravity anomaly,
            !>    no real topography)
            GrvKernel(j,i)       =
     &        1e5*((grv_uni_const*Emass)/r(nr+1)**3)*real(j-1, kind=dp)
     &        *(visc0/(bgrho*r(nr+1)*grav(r(nr+1))))*u(3)
            !> surface topography kernel
            SurKernel(j,i)=
            !  &  -(visc0/((rho0(nr+1)-rho_w)*r(nr+1)*grav(r(nr+1))))
     &     -(visc0/((bgrho-rho_w)*r(nr+1)*grav(r(nr+1))))
     &         *(u(1)+(rho_w/bgrho)*u(3))
            !> lateral velocity kernel at the surface 
            VelKernel(j,i)=
     &         u(2)*sqrt(real(j,kind=qp)*real(j+1.0, kind=qp))
            !> CMB topography kernel
            CMBKernel(j,i)=
     &       - (visc0/((rho0(1)-rho_c)*r(1)*grav(r(1))))
     &         *(u(5)+(rho_c/bgrho)*u(7))
            !> impedence/Admittance kernel Kernel
            ADMKernel(j,i) = 
     &        -(GrvKernel(j,i)/abs((SurKernel(j,i)*1e-3)))
         enddo
         !> predefine admittance kernel at CMB and surface to avoid
         !>    strange results
         ADMKernel(:,nr+1) = 0.0
         ADMKernel(:,   1) = 0.0
         if (j==1) then
            GeoKernel(j,:) = 0.0
         endif
         if (kernel_output_flg ==1) then
            call kernel_out(GeoKernel(j,:)  , 1, j)
            call kernel_out(VelKernel(j,:)  , 2, j)
            call kernel_out(CMBKernel(j,:)  , 3, j)
            call kernel_out(SurKernel(j,:)  , 4, j)
            call kernel_out(GrvKernel(j,:)  , 5, j)
            call kernel_out(ADMKernel(j,:)  , 6, j)
         end if
      endif
      enddo

      call mpi_synch_kernel(GeoKernel)
      call mpi_synch_kernel(VelKernel)
      call mpi_synch_kernel(CMBKernel)
      call mpi_synch_kernel(SurKernel)
      call mpi_synch_kernel(GrvKernel)
      call mpi_synch_kernel(ADMKernel)
      call fldcalc(GeoKernel, GrvKernel, VelKernel,
     &   SurKernel , CMBKernel)
      
      if ( inv_flg==0 .and. myproc==root_proc) then
         call summary
      endif
      end subroutine prop_matrix

*dk eigen
      subroutine eigen(a_in, n, wr_out, wi_out)
      use constants
      implicit none
      integer, parameter :: nb=64
      integer, parameter :: nmax=10
      integer, parameter :: lda=nmax
      integer, parameter :: ldvr=nmax
      integer, parameter :: lwork=(2+nb)*nmax
      integer  i, info, j, lwkopt, n
      real(qp) a_in(n,n)
      real(qp) wr_out(n), wi_out(n)
      real(dp) a(lda, nmax), dummy(1,1), vr(ldvr, nmax),
     &               wi(nmax), wr(nmax), work(lwork)

      external dgeev

      do i=1,n
         do j=1,n
            a(i,j)= real(a_in(i,j), kind=dp)
         enddo
      enddo
c..   Use Lapack
      call dgeev('No left vectors','No right Vectors',n,a,LDA,wr,wi,
     &      dummy, 1, vr, ldvr, work, lwork, info)
      if (info==0) then
         do i=1,n
            wr_out(i) = real(wr(i), kind=qp)
            wi_out(i) = real(wi(i), kind=qp)
         enddo
      else
         write(*,'(a)') 'Error in dgeev!'
         stop
      endif
      end subroutine eigen



*dk aconstruct
c.. Constructs the A matrix for a given layer
c..   and a given degree l of SPH
      subroutine aconstruct(a, ind, l)
      use constants
      use util, only : r, rho0, rvsc
      implicit none
      integer ind, i, j
      real(qp) kappa
      real(qp) grav
      real(qp) fac_visc, l, fac_rho0, fac_grav
      real(qp) a(6,6)
c..   Assumption: Viscosity of each layer is the average
c..      of the viscosities above and at the bottom
! (WARNING) SHOULD THIS BE AN AVERAGE OR THE VALUE???
      fac_visc = (rvsc(ind))! + rvsc(ind+1) )/2
      fac_rho0 = (rho0(ind))! + rho0(ind+1) )/2
	  fac_grav = (grav(r(ind)))
      do i=1,6,1
         do j=1,6,1
            a(i,j)=0.0
         enddo
      enddo
c..  See Corrieu, Thoraval & Ricard 1994 
c..   Line : 1    As a result of Continuity Equation
      a(1,1)=-2-kappa(ind)
      a(1,2)=l*(l+1)
c..   Line : 2    As a result of \theta component of 
c..               the continuity equation
      a(2,1)=-1
      a(2,2)=1
      a(2,4)=1/fac_visc
c..   Line : 3    \r component of the Stokes Equation
c..               combined with constitutive equation for 
c..               components \tau_{\theta \theta} and  \tau_{\phi \phi}
c..				  Include expression for advection of 660 km (Le Stunff 
c..				  & Ricard, 1997)
c..				  N.B. Questionable whether need reformulation based on
c..				  compressibility (see discussion in Corrieu et al. 1994)?
c..				  SEEMS TO WORK IF FAC_AD MUCH SMALLER THAN 3e10
	  if ((r(ind).ge.5700e3).and.(r(ind).le.5960e3)) then
	      a(3,1)=(4*(3+kappa(ind))*fac_visc)-
     &       ((fac_ad*fac_grav*(r(ind)**2))/visc0)
	  else
	      a(3,1)=(4*(3+kappa(ind))*fac_visc)
	  end if
      a(3,2)=-6*l*(l+1)*fac_visc
      a(3,3)=1
      a(3,4)=l*(l+1)
      a(3,6)=-fac_rho0/bgrho
c..   Line : 4    \theta component of the Stokes Equation
c..               combined with constitutive equation for 
c..               components \tau_{\theta \theta} and  \tau_{\phi \phi}
      a(4,1)=-2*(3+kappa(ind))*fac_visc
      a(4,2)=2*(2*l*(l+1)-1)*fac_visc
      a(4,3)=-1
      a(4,4)=-2
      a(4,5)=-fac_rho0/bgrho
c..   Line : 5 & 6   Simple laws of gravity and potential and the pertubation
c..                  In the potential field
      a(5,5)=1
      a(5,6)=1
	  if ((r(ind).ge.5700e3).and.(r(ind).le.5960e3)) then
     	 a(6,1)=(3*gearth*fac_ad*bgrho*(r(ind)**3))/
     &       (visc0*rhoearth*Rearth)
	  end if
      a(6,5)=l*(l+1)
      end subroutine aconstruct


*dk derivation
c..
      subroutine derivation(u_in, r_in,  u_der)
      use constants
      implicit none
      integer i
      real(qp) u_in(nr+1), r_in(nr+1)
      real(qp) u_der(nr+1)
      do i=2,nr,1
         u_der(i) = (u_in(i+1)-u_in(i-1))/(r_in(i+1)-r_in(i-1))
      enddo
      u_der(1) = 0.0
      u_der(nr+1) = 0.0
      end subroutine derivation

*dk kappa_calculation
c     Computes the kappa factor for the A matrix at each layer
      function kappa(ind) result(kappa_val)
      use constants
      use util, only : r, rho0
      implicit none
      integer ind
      real(qp) kappa_val
      if (kappaindex==1) then
         if (ind == 1 ) then
            kappa_val = 0.0
         else 
            kappa_val = (r(ind)/(rho0(ind))) * ((rho0(ind)-rho0(ind-1))/
     &       (r(ind)-r(ind-1)))
         endif
      else 
         kappa_val = 0.0
      endif
      return
      end function

