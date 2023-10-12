module sph_tools
   use constants

   implicit none


   public :: init_sph_tls, theta, phi, facto, leg 
   
   real(dp), protected, dimension(grid ), save                             :: theta
   real(dp), protected, dimension(grid2), save                             :: phi
   real(dp), protected, dimension(0:2*l_max_calc), save                    :: facto
   real(dp), protected, dimension(l_min:l_max_calc, 0:l_max_calc, grid), save  :: leg

contains

   subroutine init_sph_tls
   use SHTOOLS

   implicit none

   integer :: i, l, m, k
   real(dp) p((l_max_calc+1)*(l_max_calc+2)/2)
   real(dp) z_arg 

   leg(:,:,:) = 0.0
   do i=1,grid
      !Create colatitudes in radian based on the defined grid
      theta(i)=(90.0-(lat1+(grid-i)*Vspac))*deg2rad
      ! use SHTOOLS library to compute the associated legendre polynomials
      ! for each cos(theta(i))
      z_arg = cos(theta(i))
      call PLegendreA(p, lmax=l_max_calc, z=z_arg, csphase=1)
      ! copy legendre polynomals to our big collection of leg for all colatitudes
      k=1
      do l=0,l_max_calc
         do m=0,l
            if (l >= l_min) leg(l,m,i) = p(k)
            k=k+1
         enddo
      enddo
   enddo

   do i=1,grid2
       phi(i)=(lon1+(i-1)*Hspac)*deg2rad
   enddo

   ! all square root factorials and square roots
   facto(:)=1.0
   do l=1,2*l_max_calc
      facto(l)=l*facto(l-1)
   enddo
   facto(:)=sqrt(facto(:))

   end subroutine init_sph_tls 
end module sph_tools
