module util 
   use constants
   implicit none
contains
      !*dk initialize
      !c..   Using read_in_radial reads in the radial terms
      subroutine initialize
      use constants
      implicit none
      common /path/ output, input, name_prof
      integer(i4b) :: ir, l, m, dump
      integer(i4b) :: degmod
      real(dp), dimension(l_max_calc+1, 2*l_max_calc+1, nr+1) ::  test_anom
      real(dp) :: test1, test2
      character*100 line_temp
      character*120 output, input, name_prof
      open(001, file= trim(input), status='old', action='read')
      read(001, '(A)') line_temp
      read(line_temp(8:11) ,* ) degmod
      read(line_temp(16:19),*) dump   
      if ((dump .ne. nr) ) then
         write(*,'(A)') &
     & 'EROR: # of layers in input file and size.h inconsistent'
         stop
      endif
      read(001,*) line_temp
      read(001,*) line_temp
      call read_in_radial(r,001)
      read(001,*) line_temp
      call read_in_radial(rho0,001)
      read(001, '(A)') line_temp
      ! initialize all coeffs with zero
      anom(:,:,:) = 0.0
      do ir = nr+1,1,-1
         read(001,'(A)') line_temp
         do l=0,l_max 
            do m=0,l
               read(001,*) test1, test2
               if (l <= l_max_calc) then
               test_anom(l+1,l_max_calc+1+m,ir) = test2
               test_anom(l+1,l_max_calc+1-m,ir) = test1
               endif
            enddo
         enddo
         if ((r(nr+1)-r(ir))<(top_depth*1000) .or. &
     &       (  r(ir)-r( 1))<(bot_depth*1000)   )&
     &         test_anom(:,:,ir)=0.0
      enddo
      anom = test_anom
      close(001)
      anom(:,:,:) = anom(:,:,:)*sqt
!     TO DO : HERE MAY COME OTHER READINGS
      end subroutine


      !dk set_simple_visc
      subroutine set_simple_visc(asth_visc, asth_thick, ir_lith)
      use constants
      implicit none
      integer(i4b)      :: ir_lith, ir_asth
      real(qp)          :: asth_visc, asth_thick
      rvsc(:) = visc_mantle
      rvsc(ir_lith:nr+1) = visc_lith
      ir_asth = minloc(abs(r(:) - (r(ir_lith) - asth_thick)),1)
      rvsc(ir_asth:ir_lith-1) = asth_visc 
      end subroutine set_simple_visc

      !dk linspace
      subroutine linspace(x, x_start, x_end, x_len)
      use constants
      implicit none
      integer :: x_len, i
      real(qp) :: x_start, x_end, dx
      real(qp), dimension(x_len), intent(out) :: x
      dx = (x_end - x_start)/(x_len-1)
      x(1:x_len) = [(x_start + ((i-1)*dx),i=1,x_len)]
      end subroutine linspace

      !dk progress
      subroutine progress(j)
      implicit none
      integer(kind=4)::j,k
      character(len=37)::bar="???% |                              |"
      write(unit=bar(1:3),fmt="(i3)") 10*j/3
      do k=1, j
        bar(6+k:6+k)="*"
      enddo
      ! print the progress bar.
      write(unit=6,fmt="(a1,a1,a37)") '+',char(13), bar
      return
      end subroutine progress



      !*dk read_in_radial
      !c..   Read in radial profile u
      subroutine read_in_radial(u, in_unit)
      use constants
      implicit none                       

      integer in_unit, i
   
      real(qp) u(nr+1)
      ! Be carfule of the indexing, increament
      do i=nr+1,1,-1
         read(in_unit, *) u(i)
      enddo
      end subroutine read_in_radial


      !*dk vector_out
      !c..   writes a ni vector on the iunit
      subroutine vector_out(A,ni, iunit)
      use constants
      implicit none

      integer j
      integer ni, iunit
      real(qp) A(ni)
      write(iunit,'(A)', advance='yes') ' '
      do j=1,ni
         write(iunit, 10, advance='no') A(j)
      enddo
      write(iunit,'(A)', advance='yes') ' '

 10   format(f9.5,2x)
      end subroutine vector_out

      !*dk matrix_out
      !c..   writes a ni*nk matrix on the iunit
      subroutine matrix_out(A,ni,nk, iunit)
      use constants
      implicit none

      integer j,k
      integer nk, ni, iunit
      real(qp) A(ni, nk)
      do j=1,ni
         do k=1,nk
            write(iunit, 10, advance='no') A(j,k)
         enddo
         write(iunit,'(A)', advance='yes') ' '
      enddo
 10   format(e10.3,1x)
      end subroutine matrix_out

      !*dk kernel_out
      subroutine kernel_out( u, io , l)
      implicit none
      common /path/ output, input, name_prof
      integer l, i, io
      real(dp) u(nr+1)
      character*20 char1, char2
      character*120 output, input, name_prof
      if (io .eq. 1) then
         char1 = 'geoid'
      elseif (io .eq. 2) then
         char1 = 'surfvel'
      elseif (io .eq. 3) then
         char1 = 'cmbtopo'
      elseif (io .eq. 4) then
         char1 = 'surftopo'
      elseif (io .eq. 5) then
         char1 = 'gravity'
      elseif (io .eq. 6) then
         char1 = 'admittance'
      else
         write(*,*) 'No valid io for kernel_out!'
         stop
      endif
      call system('mkdir -p '//trim(output)//'KERNELS/')
      write(char2, '(I3.3)') l
      open(333, file=trim(output)//'KERNELS/'//trim(char1)&
     & //trim(char2), action = 'write', status='unknown')
      do i=1,nr+1,1
         write(333, '(f17.2,4X, es16.7e3)') r(i)/1000, u(i)
      enddo
      close(333)
      end subroutine kernel_out

      subroutine sphcout_3d(coef, char1)
      common /path/ output, input, name_prof
      integer(i4b) ir, l, m
      real(dp) coef (0:l_max_calc, -l_max_calc:+l_max_calc,1:nr+1)
      character*15 char1
      character*120 output, input, name_prof
      open(002,file=trim(output)//'/'//trim(char1),action='write')
      write(002,'(A,I4,A,I4,A)') '# LMAX=',l_max_calc, ', NR',nr,' 1'
      write(002,'(A,A)') '# Temperature: INPUT : ', trim(input)
      write(002,'(A,A)') '# Radius'
      do ir= nr+1,1,-1
         write(002,'(f15.2)') r(ir)
      enddo
      write(002,'(A,A)') '# Averages'
      do ir= nr+1,1,-1
         write(002,'(f10.2)') coef(0,0,ir)
      enddo
      write(002,'(A,A)') '# Spherical Harmonics'
      do ir= nr+1,1,-1
         write(002,'(A,I4.1)') '# Radial Index:', ir 
          do l=0,l_max_calc
             do m=0,l
                if (m==0) then
                   write(002,'(e15.6, e15.6)')  coef(l,m,ir), 0.0
               else
                  write(002,'(e15.6, e15.6)') coef(l,-m,ir), coef(l,+m,ir)
                endif
             enddo
          enddo
      enddo
      close(002)
      end subroutine sphcout_3d

      subroutine sphcout(coef,char1)
      implicit none
      common /path/ output, input, name_prof
      integer(i4b) l, m
      real(dp) coef (l_min:l_max_calc, -l_max_calc:+l_max_calc)
      character(len=*) :: char1
      character*120 output, input, name_prof
      open(444, file=trim(output)//trim(char1),&
     & action='write', status='unknown')
!     The header part
      write(444,'(A)') '... HEADER'
      write(444,'(A,A)') '... INPUT : ', trim(input)
      write(444,'(A,A)') '... ', trim(output)
      write(444,'(A)') 'l ,m , Cnm(m>=0) , Snm(m<0)'
      do l=l_min,l_max_calc
         do m=0,l
            if (m==0) then
                write(444,'(A,I8,I8,e15.6, e15.6)') 'SPHCO',l, m, coef(l,m), 0.
            else
                write(444,'(A,I8,I8,e15.6, e15.6)')  'SPHCO',l, m, coef(l,-m), coef(l,+m)
            endif
         enddo
      enddo
      write(444,'(A)') '... END OF FILE!'
      close(444)
      end subroutine sphcout

      subroutine inv_out_rslt(p, fret, iter)
      use constants , only : np
      use nrtype, only : sp, i4b
      implicit none
      common /path/ output, input, name_prof
      integer(i4b) :: out_unt, ip, ir
      integer(i4b),  intent(in)                :: iter
      real(sp),      intent(in)                :: fret
      real(sp),      intent(in), dimension(np) :: p
      character*120 output, input, name_prof
      character*3 step_char

      write(step_char,'(I3.3)') iter
      out_unt = get_file_unit(999)
      open(unit= out_unt,&
     & file=trim(output)//'/powell_prm_iter_'//step_char)
      write(out_unt, '(I3.1, E10.3)') iter, fret
      write(out_unt, '(A)')  '-----------------'
      do ip=1,np
         write(out_unt, '(E10.3)') p(ip)
      enddo
      close(out_unt)
      open(unit= out_unt,&
     & file=trim(output)//'/powell_visc_iter_'//step_char)
      write(out_unt, '(I3.1, E10.3)') iter, fret
      write(out_unt, '(A)')  '-----------------'
      do ir=nr+1,1, -1
         write(out_unt, '(F10.2, E10.3)') r(ir), rvsc(ir)
      enddo
      close(out_unt)
      end subroutine inv_out_rslt

      subroutine fldout(field, char1)
      use constants
      use sph_tools, only : theta, phi
      implicit none
      common /path/ output, input, name_prof
      real(dp) field(grid, grid2)
      integer(i4b) :: i, j, fldoutid
      character(len=*) :: char1
      character*120 output, input, name_prof
      open(645, file=trim(output)//trim(char1), status = "unknown", action = 'write')
      do i=1,grid
         do j=1,grid2
            write(645,'(f7.2,3X,f7.2,3X,f15.7)') phi(j)*rad2deg,  (pi/2-theta(i))*rad2deg, field(i,j)
         enddo
      enddo
      close(645)
      end subroutine fldout

      !*dk summary 
      subroutine summary
      use constants
      implicit none
      common /path/ output, input, name_prof
      common /ResCorr/ GeoCorr, GravCorr, DyntCorr,CMBCorr, surmax, geomax
      integer(i4b) :: i,l,m
      character*120 output, input, name_prof
      real(qp) kappa
      real(qp) grav
      real(dp) GeoCorr, GravCorr, DyntCorr, CMBCorr, surmax, geomax
      open(111, file=trim(output)//'Summary', action = 'write', status = 'unknown' )
      write(111, '(5X,A,A)') 'Summary: Input File - ', trim(input)
      write(111, '(15X, A)') 'Saved to: ',trim(output)
      do i = 1,62
         write(111,'(A)', advance = 'no') '*'
      enddo
      write(111,*) 
      write(111, 22)nr*2, nr+1, l_min, l_max_calc, top_depth, bot_depth,&
     & bdry, rho_c, lat1, lat2, lon1, lon2, GeoCorr, geomax, surmax,&
     & GravCorr
 22   format &
     &(' MT    =',I10  ,'| NR+1    =',I10  ,'| LMIN    = ',I10  ,/,&
     & ' LMAX  =',I10  ,'| TOP_DPTH=',f10.2,'| BOT_DPTH= ',f10.2,/,&
     & ' bdry  =',I10  ,'| RHO_C   =',e10.2,'| LAT1    = ',f10.2,/,&
     & ' LAT2  =',f10.2,'| LON1    =',f10.2,'| LON2    = ',f10.2,/,&
     & 'GeoCorr=',f10.2,'| MaxGeoid=',f10.2,'| MaxSurTo= ',f10.2,/,&
     & 'GrvCorr=',f10.2)
      do i = 1,62
         write(111,'(A)', advance = 'no') '*'
      enddo
      write(111,*) 
      write(111, '(3x,3A,A,A)')'Radius ',' Viscosity ' &
     &,' Density ',' CHI ',' Gravity ' 
      do i = nr+1, 1, -1        
         write(111,'(f15.0, e10.3, f8.2, e10.1, f7.2)') &
     &            r(i),rvsc(i)*visc0, rho0(i), kappa(i), grav(r(i))
      enddo
      close(111)
      open(111,file=trim(output)//'Density', action = 'write',&
     &     status = 'unknown')
      
      write(111,'(A)') ' l  m     Radius   Cnm   Snm'
      do i=nr+1,1,-1
         do l=0,l_max_calc
            do m=0,l
               if (m==0) then
                  write(111,'(I2.2,3X,I2.2,3X,I3.3,3X,f11.3,3X,E12.4,3X,E12.4)') &
     &      l,m, i, r(i),anom(l,m,i)/sqt, 0.00          
               else                                    
                  write(111,'(I2.2,3X,I2.2,3X,I3.3,3X,f11.3,3X,E12.4,3X,E12.4)')&
     &      l,m,i, r(i),anom(l,m,i)/sqt,anom(l,-m,i)/sqt
               endif
            enddo
         enddo
      enddo

      close(111)
      end subroutine summary

      !*dk viscinit
      subroutine viscinit
      use constants
      implicit none
      common /path/ output, input, name_prof
      integer(i4b) ::  i
      character*120 input, output, name_prof
      open(183, file=trim(name_prof),action='read', status='old')
      do i=nr+1,1,-1
         read(183,*) rvsc(i)
      enddo
      close(183)
      end subroutine viscinit

      function get_file_unit(lu_max) result(out_unit)
      use constants, only : i4b
      implicit none
      ! get_file_unit returns a unit number that is not in use
      integer(i4b), intent(in) :: lu_max
      integer(i4b) :: out_unit
      integer(i4b) :: lu, m, iostat
      logical ::   opened
      m = lu_max  ;  if (m < 1) m = 97
      do lu = m,1,-1
         inquire (unit=lu, opened=opened, iostat=iostat)
         if (iostat.ne.0) cycle
         if (.not.opened) exit
      end do
      out_unit = lu
      return
      end function get_file_unit
   ! ------------------------------------------------------------
   !               convolution of two functions x and h
   ! ------------------------------------------------------------
   function convolve(x, h)
      use constants, only : qp, i4b
      implicit none
      
      !x is the signal array
      !h is the noise/impulse array
      real(qp), dimension(:), allocatable :: convolve, y
      real(qp), dimension(:) :: x, h
      integer(i4b) :: kernelsize, datasize
      integer(i4b) :: i,j,k, kk
      
      datasize       = size(x)
      kernelsize     = size(h)
      allocate(y        (datasize+kernelsize))
      allocate(convolve (datasize))
      y(:) = 0.0
      !first part
      do i= 1, datasize
         kk = 0
         do k= i, 1, -1
            if (kernelsize-kk .gt. 0)  then
               y(i) = y(i) + x(k)*h(kernelsize-kk)
            endif
            kk = kk + 1
         enddo
      enddo
      !last part
      do i= 1, kernelsize 
         kk = 1
         do k= datasize-kernelsize+i+1, datasize
            y(datasize+i) = y(datasize+i) + x(k)*h(kk)
            kk = kk +1
         enddo
      enddo
      y(:) = y(:)/maxval(y)*maxval(x)
      convolve(:) = y(int(kernelsize/2)+1:datasize+int(kernelsize/2))
    end function convolve
   ! --------------------------------------------------------------------
   !     gaussian function
   ! --------------------------------------------------------------------
   function gaussian(sigma, nx)
      use constants, only : i4b, dp, qp, pi
      implicit none
      integer(i4b)               :: nx, ix
      real(qp)                   :: y_bound
      real(qp)                   :: sigma
      real(qp), dimension(nx+1)  :: gaussian
      y_bound = real(nx, kind=qp)
      do ix = 1, nx
         gaussian(ix) = -1*y_bound + (ix-1)*((2*y_bound)/(nx-1.0))
      enddo
      gaussian = 1/sqrt(2*pi*sigma**2)*exp(-gaussian**2/(2*sigma**2))
      return
   end function gaussian
   subroutine meanvariance(in_arr, n, stddev)
   ! --------------------------------------------------------------------
   ! meanvariance:
   !    this routine reads in an array and computes the mean, variance
   ! and standard deviation of the tomographic_models stored in the array.
   ! --------------------------------------------------------------------
      use constants, only : qp, dp
      implicit  none
      integer(i4b)             :: n          ! actual array size
      real(qp), dimension(1:n) :: in_arr     ! input array
      real(qp)                 :: mean, variance
      real(dp)                 :: stddev        ! results
      integer(i4b)             :: i                   ! running index
      mean = 0.0                           ! compute mean
      do i = 1, n
         mean = mean + in_arr(i)
      end do
      mean = mean / n
      variance = 0.0                       ! compute variance
      do i = 1, n
         variance = variance + (in_arr(i) - mean)**2
      end do
      variance = variance / (n - 1)
      stddev   = real(sqrt(variance), kind=dp)            ! compute standard deviation
   end subroutine meanvariance
end module util

