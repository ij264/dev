module constants 
   implicit none
   include "constants.h"
   real(qp),       dimension(nr+1)   :: r, rho0
   real(qp),       dimension(nr+1)   :: rvsc 
   real(dp),       dimension(0:l_max_calc,-l_max_calc:+l_max_calc,1:nr+1) ::  anom
   data prms_bounds_dpth /3000e3, 1100e3, 600e3, 0e3/ 
   data prms_begin       / 5e6, 5e5,  1e5/ 
contains
   subroutine init_const 
   implicit none
   integer(i4b) :: ir

   prms_bounds_rshl(:) = Rearth - prms_bounds_dpth(:)

   do ir=1,np+1
      prms_bounds_indx(ir) = int(minloc(abs(r-prms_bounds_rshl(ir)),dim=1), kind=i4b)
   enddo
   end subroutine init_const 
end  module constants
