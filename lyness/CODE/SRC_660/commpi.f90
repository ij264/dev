module commpi
   use mpi
   implicit none
!    include "mpif.h"
   integer :: myproc, nproc, ierror
   integer :: dp_data_type=MPI_REAL8 
   integer :: root_proc = 0               
contains

   !dk pinit
   subroutine pinit
   ! This routine initializes the MPI processes.
      implicit none
      
      !Enroll this process in MPI.
      call mpi_init( ierror )
   
      call mpi_comm_rank( mpi_comm_world, myproc, ierror )
      call mpi_comm_size( mpi_comm_world, nproc, ierror )
      call mpi_type_create_f90_real(15, 307, dp_data_type, ierror )
   end subroutine pinit

   !dk pinit
   subroutine wait4others
      implicit none
      call mpi_barrier( mpi_comm_world, ierror )
   end subroutine wait4others
   
   !*dk pexit
   subroutine pexit
  
      call mpi_barrier( mpi_comm_world, ierror )
      !if (myproc .neqv. root_proc)
      !call mpi_abort(ierror)
      call mpi_finalize( ierror )
   
   end subroutine pexit


   subroutine mpi_synch_coef(coef)
      use constants
      implicit none
      integer ierr
      real(dp) coef        (l_min:l_max_calc, -l_max_calc:+l_max_calc) 
      real(dp) dummy_send  (l_min:l_max_calc, -l_max_calc:+l_max_calc) 
      real(dp) dummy_rec   (l_min:l_max_calc, -l_max_calc:+l_max_calc) 
      
      dummy_send = coef
      call mpi_allreduce(dummy_send, dummy_rec, (l_max_calc+1)*(2*l_max_calc+1), mpi_real8, mpi_sum, mpi_comm_world, ierror)
      coef = dummy_rec
   end subroutine mpi_synch_coef



   subroutine mpi_synch_fld(fld)
      use constants 
      implicit none
      integer ierr
      real(dp) fld         (grid, grid2) 
      real(dp) dummy_send  (grid, grid2) 
      real(dp) dummy_rec   (grid, grid2) 
      
      dummy_send = fld
      call mpi_allreduce(dummy_send, dummy_rec, grid*grid2, mpi_real8, mpi_sum, mpi_comm_world, ierror)
      fld = dummy_rec
   end subroutine mpi_synch_fld


   !*dk mpi_synch_kernel
   subroutine mpi_synch_kernel(all_kernels)
      use constants   
      implicit none
   
      integer j, ierr
      real(dp) all_kernels(l_min:l_max_calc, nr+1) 
      real(dp) dummy_send (l_min:l_max_calc, nr+1) 
      real(dp) dummy_rec  (l_min:l_max_calc, nr+1) 
      
      dummy_send = all_kernels
      call mpi_allreduce(dummy_send, dummy_rec, (l_max_calc-l_min+1)*(nr+1), mpi_real8, mpi_sum, mpi_comm_world, ierror)
      all_kernels = dummy_rec
   end subroutine mpi_synch_kernel
end module commpi
