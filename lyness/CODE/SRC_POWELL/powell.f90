   subroutine powell(p,xi,ftol,iter,fret)
   use commpi
   use constants, only : np, inv_out_skp
   use nrtype; use nrutil, only : assert_eq,nrerror
   use util, only : inv_out_rslt 
   use nr, only : linmin
   implicit none
   real(sp), dimension(np), intent(inout) :: p
   real(sp), dimension(np,np), intent(inout) :: xi
   integer(i4b), intent(out) :: iter
   integer(i4b) :: ip
   real(sp), intent(in) :: ftol
   real(sp), intent(out) :: fret
   character*3 step_char
   interface
      function misfit_func(p)
      use nrtype
      use constants, only : np
      implicit none
      real(sp), dimension(np), intent(in) :: p
      real(sp) :: misfit_func
      end function misfit_func
   end interface
   integer(i4b), parameter :: itmax=200
   real(sp), parameter :: tiny=1.0e-25_sp
   integer(i4b) :: i,ibig,n
   real(sp) :: del,fp,fptt,t
   real(sp), dimension(size(p)) :: pt,ptt,xit
   n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
   fret=misfit_func(p)
   pt(:)=p(:)
   if (myproc==root_proc) print *, 'Start of Powell Minimization.'
   iter=0
   do
      if (myproc==root_proc) then 
      write(*,'(A,I3.1,A,e12.5)') "Iter: ",iter, " Current MINIMUM: ", fret
      if (mod(iter, inv_out_skp+1)==0) then
         call inv_out_rslt(p, fret, iter)
      endif
      endif
      iter=iter+1
      fp=fret
      ibig=0
      del=0.0
      do i=1,n
         xit(:)=xi(:,i)
         fptt=fret
         call linmin(p,xit,fret)
         if (fptt-fret > del) then
            del=fptt-fret
            ibig=i
         end if
      end do
      if (2.0_sp*(fp-fret) <= ftol*(abs(fp)+abs(fret))+tiny) return
      if (iter == itmax) call &
         nrerror('powell exceeding maximum iterations')
      ptt(:)=2.0_sp*p(:)-pt(:)
      xit(:)=p(:)-pt(:)
      pt(:)=p(:)
      fptt=misfit_func(ptt)
      if (fptt >= fp) cycle
      t=2.0_sp*(fp-2.0_sp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if (t >= 0.0) cycle
      call linmin(p,xit,fret)
      xi(:,ibig)=xi(:,n)
      xi(:,n)=xit(:)
   end do
   end subroutine powell
