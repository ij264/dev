module f1dim_mod
   use nrtype
   integer(i4b) :: ncom
   real(sp), dimension(:), pointer :: pcom,xicom
contains
!bl
   function f1dim(x)
   implicit none
   real(sp), intent(in) :: x
   real(sp) :: f1dim
   interface
      function misfit_func(x)
      use nrtype
      use constants, only : np
      real(sp), dimension(np), intent(in) :: x
      real(sp) :: misfit_func
      end function misfit_func
   end interface
   real(sp), dimension(:), allocatable :: xt
   allocate(xt(ncom))
   xt(:)=pcom(:)+x*xicom(:)
   f1dim=misfit_func(xt)
   deallocate(xt)
   end function f1dim
end module f1dim_mod

   subroutine linmin(p,xi,fret)
   use nrtype; use nrutil, only : assert_eq
   use nr, only : mnbrak,brent
   use f1dim_mod
   implicit none
   real(sp), intent(out) :: fret
   real(sp), dimension(:), target, intent(inout) :: p,xi
   real(sp), parameter :: tol=1.0e-4_sp
   real(sp) :: ax,bx,fa,fb,fx,xmin,xx
   ncom=assert_eq(size(p),size(xi),'linmin')
   pcom=>p
   xicom=>xi
   ax=0.0
   xx=1.0
   call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
   fret=brent(ax,xx,bx,f1dim,tol,xmin)
   xi=xmin*xi
   p=p+xi
   end subroutine linmin
