      program driver
      use constants
      use ref_flds, only : ref_flds_init
      use sph_tools, only : init_sph_tls, theta, phi
      use util, only : r, rho0, rvsc, initialize, viscinit
      use commpi
      implicit none

      common /path/     output, input, name_prof

      interface 
         subroutine inviscpll
         end subroutine inviscpll
      end interface

      integer N, sufi, suff, step
      real(dp) time_1, time_2
      real(dp) start, finish
      character*200 output, input, name_prof


      N = IARGC()

      call cpu_time(time_1)
      if (N .ne. 2 ) then
         write(*,'(A)') 'Number of given arguments should be 3!'
         write(*,'(I1.1, A)') N, ' is given. Try the following order:'
         write(*,'(A)') './GEOID [OUTPUT] [VISC_PROFILE]'
         stop
      endif

      call getarg(1,output)      ! ouput is the path to output data
      call getarg(2,name_prof)   ! viscosity profile
      input='inputs.sph'
      
      call cpu_time(start)
      
      call system('mkdir -p '//trim(output))

      ! initialize the parallel process
      call pinit

      !> reading the radius
      call initialize
      !> initialise the reference fields
      call ref_flds_init
      !> initialise the legendre polynomials
      call init_sph_tls 
      !> initialize the necessary parameters and constants, in module consts
      call init_const
      !> initialize viscosity profile
      call viscinit
      !> compute kernels using propagator matrix tech.
      call prop_matrix
      !> write out the computed correlation
      if (CorrFlag == 1 .and. myproc==root_proc) then
          write(*,'(A)') '-------------------------------------'
          write(*,'(A, I2.1)') ' Processor reporting:', myproc
          write(*,'(A)') '-------------------------------------'
          call cpu_time(finish)
      endif

      call cpu_time(time_2)
      if (myproc==root_proc) write(6,'(A,I2.2,A,I2.2,A,I2.2,A)')
     & '  Elapsed Time:                 ',
     & int((time_2-time_1)/3600), ':',
     &   int((time_2 - time_1)/60),''':',mod(int(time_2-time_1),60),'"'

      call pexit
      end program driver

*dk subroutine prop_matrix
      subroutine prop_matrix

      use commpi
      use constants
      use util, only : r, rho0, rvsc, kernel_out, summary
      implicit none
      common /path/ output, input, name_prof
      integer i, j, k, m, n, ierr, iiir
      integer ipiv(8), kki, kkj
      character*200 output, input, name_prof
      real(dp) GeoKernel(l_min:l_max_calc,nr+1)
      real(dp) VelKernel(l_min:l_max_calc,nr+1)
      real(dp) SurKernel(l_min:l_max_calc,nr+1)
      real(dp) GrvKernel(l_min:l_max_calc,nr+1)
      real(dp) ADMKernel(l_min:l_max_calc,nr+1)
      real(dp) CMBKernel(l_min:l_max_calc,nr+1)
      real(qp) kappa
      real(qp) a(6,6), prop(6,6,1:nr+2)
	  real(qp) afull(6,6,1:nr+2)
      real(qp) eigen_i(6), eigen_r(6)
      real(qp) c(8,8), b(8), u(8)
      real(qp) dump, fac
      real(qp) grav
      real(qp) rcond
      real(qp) geopref
      real(qp) surfgrav
      real(qp) DTpref
      real(qp) rhom ! mantle density
      real(qp) rhow ! water density

      rhom = 3380
      rhow = 1030

      GeoKernel(:,:) = 0.0
      VelKernel(:,:) = 0.0
      CMBKernel(:,:) = 0.0
      SurKernel(:,:) = 0.0
      GrvKernel(:,:) = 0.0
      ADMKernel(:,:) = 0.0

      !> compute the gravity at the surface
      surfgrav = grav(r(nr+1))

      !> constructing the matrix
      do j= l_min, l_max_calc, 1
      if (mod(j-l_min,nproc) == myproc) then
         !> initialisation of a 6*6 matrix with zero values
         call matunitinit(prop(:,:,nr+2), 6,6)
         do i=nr+1,2,-1
            call aconstruct(a, i, real(j, kind=qp))
            call eigen(a, 6, eigen_r, eigen_i)
            call propmat(a, prop(:,:,i), eigen_r,eigen_i,r(i)/r(i-1))
            prop(:,:,i) = matmul(prop(:,:,i+1), prop(:,:,i))
         enddo
         !> prop(:,:,2) represents prop{6*6}(c,a), i.e. core to surface
		 !> construct C matrix that includes boundary conditions, now 8x8.
		 !> C matrix is essentially equivalent to specific version of A matrix
		 !> that is relevant to solution from core to surface.
         call cconstruct(c, j, prop(:,:,2))
		 !> factor C matrix so it can be solved u at this point is null
         call qgeco(c,8,8, ipiv, rcond, u)

         if (rcond .lt. 0.0) then
            write(*,*) rcond
            stop
         endif

         geopref = 4*pi*grv_uni_const*Rearth/((2*j+1)*surfgrav)
         DTpref = 1/(rhom-rhow)

         do i=nr+1,1,-1
			!> construct b vector (i.e. RHS vector)
            call bconstruct(i, prop(:,:,i+1), j, u)
			!> solve the C matrix outputting u, where u is 
			!> equivalent to x in A*x = b on output, C is 
			!> equivalent to A and u is equivalent to b on input
 			call qgesl(c, 8, 8, ipiv, u, 0)

            !> geoid kernel

            GeoKernel(j,i)       =
     &        (visc0/(bgrho*r(nr+1)*surfgrav))*u(3)/geopref
            !> free air gravity anomaly kernel, (in our idealized geometry
            !>    gravity anomaly is the same as free air gravity anomaly,
            !>    no real topography)
            GrvKernel(j,i)       =
     &        1e5*((grv_uni_const*Emass)/r(nr+1)**3)*real(j-1, kind=dp)
     &        *(visc0/(bgrho*r(nr+1)*grav(r(nr+1))))*u(3)
            !> surface topography kernel
            SurKernel(j,i)=
     &  -(visc0/((rho0(nr+1)-rho_w)*r(nr+1)*grav(r(nr+1))))
C      &     -(visc0/((bgrho-rho_w)*r(nr+1)*grav(r(nr+1))))
     &         *(u(1)+(rho_w/bgrho)*u(3))/DTpref
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
      real(qp) fac_visc, l, fac_rho0
      real(qp) a(6,6)
c..   Assumption: Viscosity of each layer is the average
c..      of the viscosities above and at the bottom
! (WARNING) SHOULD THIS BE AN AVERAGE OR THE VALUE???
      fac_visc = (rvsc(ind))! + rvsc(ind+1) )/2
      fac_rho0 = (rho0(ind))! + rho0(ind+1) )/2
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
      a(3,1)=4*(3+kappa(ind))*fac_visc
      a(3,2)=-6*l*(l+1)*fac_visc
      a(3,3)=1
      a(3,4)=l*(l+1)
      a(3,6)=- fac_rho0/bgrho
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

