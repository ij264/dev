   subroutine inviscpll
   use nrtype
   use constants, only : np, prms_begin
   use util, only : inv_out_rslt
   implicit none
   real(sp), dimension(np)       :: p
   real(sp), dimension(np, np)   :: xi
   real(sp)                      :: fret
   real(sp)                      :: ftol
   integer(i4b)                  :: iter
   integer(i4b)                  :: ip
   interface 
      function misfit_func(p)
         use nrtype
         use constants, only : np
         implicit none
         real(sp), dimension(np), intent(in) :: p
         real(sp)                            :: misfit_func
      end function misfit_func
   end interface
   xi(:,:)      =  0.0 
   do ip=np, 1, -1
      xi(ip,ip) = 10.0
      p(ip)     = prms_begin(ip)
   enddo
   ftol         = 0.0001
   call powell(p, xi, ftol, iter, fret)
   call inv_out_rslt(p, fret, iter)
   end subroutine inviscpll

   function misfit_func(p)
   use constants, only : np, nr, rvsc
   use nrtype
   use commpi
   implicit none
   common /ResCorr/ GeoCorr, GravCorr, DyntCorr, CMBCorr, surmax, geomax
   common /ResStat/ flatness, obs_surmax, obs_geomax
   logical(lgt) :: err_flg
   integer(i4b) :: ii
   real(sp), dimension(np), intent(in) :: p
   real(sp)                            :: misfit_func
   real(dp)                            :: GeoCorr, GravCorr, DyntCorr, CMBCorr, surmax, geomax
   real(dp)                            :: flatness, obs_surmax, obs_geomax
   real(dp)                            :: add_misfit
   ! transparent output
   write(6,'(A)', advance='no') 'Ps to try:'
   do ii = 1, np
      write(6,'(E12.4)', ADVANCE='no') p(ii)
   enddo

   ! set viscosity based on input parameters
   call set_new_visc(p, err_flg)
   
   ! solve the equations using the new viscosity radial profile
   if (err_flg) then
      misfit_func = 1e31
   else 
      call prop_matrix
      misfit_func = ((1- (0.9*GeoCorr+0.1*CMBCorr))+ abs(geomax-obs_geomax)/100)*100
   endif
   write(6,'(A, f12.3)') ', MIS: ', misfit_func
   !misfit_func = (p(1)+2*p(2)-7)**2 + (2*p(1)+p(2)-5)**2
   end function misfit_func

   subroutine set_new_visc(p, err_flg)
   use constants
   use nrtype, only : lgt
   use util, only : gaussian, convolve
   implicit none
   integer, parameter :: sp = kind(1.0)
   integer(i4b) :: ir, kr
   logical(lgt) :: err_flg
   real(sp), dimension(np), intent(in) :: p
   real(qp), dimension(2*hlf_ny+1)     :: gaussfunc

   err_flg  = .false.
   do ir=1, np
      do kr=prms_bounds_indx(ir), prms_bounds_indx(ir+1)
         rvsc(kr) = real(p(ir), kind=qp)
         if (rvsc(kr) .lt. 1e-2) then
            err_flg = .true.
         endif
         if (rvsc(kr) .gt. 1e22) then
            err_flg = .true.
         endif
      enddo
   enddo
   rvsc(1:3) = rvsc(4)*1e-2
   rvsc(nr-1:nr+1) = rvsc(nr-2)*1e+2;
   gaussfunc = gaussian(gaus_sigma, 2*hlf_ny+1)
   rvsc = convolve(rvsc, gaussfunc)
   end subroutine set_new_visc  
