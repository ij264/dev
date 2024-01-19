*dk propmat
      subroutine propmat(a,b,u,w, rad)
      use constants
      implicit none
      integer i, k, n, m
      real(qp) rad
      real(qp) eig(6), u(6), w(6)
      real(qp) a(6,6), b(6,6), uni(6,6)
      real(qp) temp(6,6)
      do i=1,6,1
         if (w(i) .ne. 0) write(*,*) 
     & 'Warning: Non-Zero imaginary eigenvalue'
      enddo
c..   Set the repeated eigen values to zero
      eig = 0.0
      eig(1)=u(1)
      m = 1
      do i=1,6
         n = 0   
         do k=1,m
            if (abs(u(i)-eig(k)) .lt. 10E-4 ) then
                n= n+1
            endif
         enddo
         if (n .eq. 0) then
            m = m+1
            eig(m) = u(i)
         endif
      enddo
!      if (m .ne. 4) write(*,*)
!     &      'Warning: m>4, the loop is not suitable' 
!A-e(id(1))*uni,A-e(id(2))*uni),A-e(id(3))*uni)&
!     &   /(e(i)-e(id(1)))/(e(i)-e(id(2)))/(e(i)-e(id(3)))*rad**e(i)
      do i=1,6
         do k=1,6
            if (i==k) then
               uni(i,k) = 1.0
            else
               uni(i,k) = 0.0
            endif
         enddo
      enddo
      call matzeroinit(b,6,6)
      do k=1,m,1
         call matunitinit(temp,6,6)
         do i=1,m,1
            if (i .ne. k) then
               temp = matmul(temp,(a-(eig(i)*uni)))/(eig(k)-eig(i))
            endif
         enddo
         temp = temp*(rad**eig(k))
         b = b + temp
      enddo
      end subroutine

*dk matzeroinit
      subroutine matzeroinit(a,n,m)
      use constants
      implicit none
      integer n,m,i,j
      real(qp) a(n,m)
      do i=1,n,1
         do j=1,n,1
            a(i,j) = 0.0
         enddo
      enddo
      end subroutine matzeroinit

* function grav
      function grav(radius) result(grav_val)
      use constants
      use util, only : r, rho0
      implicit none
      integer(i4b) :: i
      real(qp) radius
      real(qp) grav_val
      real(qp), parameter:: fac=4*pi/3*grv_uni_const  ! 4/3*pi*G
      
      ! gravitational acc. generated by the core
      grav_val = 0.0
!      grav= fac*rho_c*(r(1)**3)/(radius**2)
!     I just put sth here, so that the grav on CMB would be consistent with PREM.csv and on surface 9.8
      grav_val = 130522366220506./(radius**2)
!      write(*,*)'here it comes',fac*rho_c*(r(1)**3)
!      stop
!      grav= fac*rho_c*(r(1)**3)/(radius**2)
      ! gravitational acc. generated by the density layers below the current position
      i = 2
      do while(r(i) .lt. radius )
          grav_val =
     & grav_val + fac*rho0(i)*((r(i)**3)-(r(i-1)**3))/(radius**2)
          i=i+1
      enddo
      ! gravitational acc. generated by the area from the
      ! nearest layer boundary to the current position
      grav_val = grav_val
     &   +fac*rho0(i)*((radius**3)-r(i-1)**3)/(radius**2)
      end function grav

* dk bconstruct
      subroutine bconstruct(ind, p, l, b)
      use constants
      use util, only : r, rho0, rvsc
      implicit none
      integer(i4b) ind, k, l
      real(qp) b(8), p(6,6)
      real(qp) const , ave_grav, ave_rad
      real(qp) grav
      const = 4*pi*grv_uni_const/(2*l+1)
      ave_rad = r(ind) !+ r(ind+1))/2
      ave_grav= grav(r(ind))! + grav(r(ind+1)))/2
      do k=1,6
         b(k)=(p(k,3)*ave_rad*ave_grav/visc0)
     &   -  p(k,6)*4*pi*ave_rad*ave_rad*grv_uni_const*bgrho/visc0
      enddo
!      b(7)=bgrho/visc0*const*((r(ind)/r(nr+1))**l)*(r(ind)**2)
!      b(8)=bgrho/visc0*const*((r(1)/r(ind))**(l-1))*(r(1)**2)
      b(7)=bgrho/visc0*const*((ave_rad/r(nr+1))**l)*(ave_rad**2)
      b(8)=bgrho/visc0*const*((r(1)/ ave_rad )**(l-1))*(r(1)**2)
      end subroutine bconstruct

*dk cconstruct
      subroutine cconstruct(C, l, p)
      ! constructs the C matrix given the p matrix
      use constants
      use util, only : r, rho0 
      implicit none
      integer(i4b) l, k, m, n
      real(qp) C(8,8), p(6,6), fac
      real(qp) grav
      call matzeroinit(C,8,8)
      fac = 4*pi*grv_uni_const/(2*real(l,kind=qp)+1)
      C(2,2)= real(bdry,kind=qp)
      C(3,1)= 1
      C(4,2)= real(1-bdry,kind=qp)
      C(5,3)= 1
      C(6,4)= 1
      C(7,1)= fac*r(nr+1)/grav(r(nr+1))*bgrho
      C(7,3)= 1
      C(7,5)= -fac*r(1)/grav(r(1))
     &         *bgrho*((r(1)/r(nr+1))**l)
      C(7,7)= C(7,5)*rho_c/bgrho
      C(8,1)= C(7,1)*((r(1)/r(nr+1))**(l+1))
      C(8,5)= -fac*r(1)/grav(r(1))*bgrho
      C(8,7)= 1-fac*r(1)/grav(r(1))*rho_c

      do k=1,6
            do m=5,8
!               n=m-3
!               if(m>6) n=n+1
               C(k,m)=-p(k,m-2)
!               C(k,m)=-p(k,n)
            enddo
      enddo
      if (bdry_bot==1) then
         do k=1,6
            C(k,6) = -p(k,2)
!            C(k,5) = -p(k,2)
!            C(k,6) = -p(k,3)
         enddo
      endif
      end subroutine

*dk solve_lgs
      subroutine solve_lgs(C, B, N, U, IPIV)
      use constants
      implicit none
      
!  Using lapack package to solve Ax=B      
!subroutine    dgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
!Parameters
!    [in] N  
!
!              N is INTEGER
!              The number of linear equations, i.e., the order of the
!              matrix A.  N >= 0.
!
!    [in] NRHS  
!
!              NRHS is INTEGER
!              The number of right hand sides, i.e., the number of columns
!              of the matrix B.  NRHS >= 0.
!
!    [in,out]   A  
!
!              A is DOUBLE PRECISION array, dimension (LDA,N)
!              On entry, the N-by-N coefficient matrix A.
!              On exit, the factors L and U from the factorization
!              A = P*L*U; the unit diagonal elements of L are not stored.
!
!    [in] LDA   
!
!              LDA is INTEGER
!              The leading dimension of the array A.  LDA >= max(1,N).
!
!    [out]   IPIV  
!
!              IPIV is INTEGER array, dimension (N)
!              The pivot indices that define the permutation matrix P;
!              row i of the matrix was interchanged with row IPIV(i).
!
!    [in,out]   B  
!
!              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!              On entry, the N-by-NRHS matrix of right hand side matrix B.
!              On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!    [in] LDB   
!
!              LDB is INTEGER
!              The leading dimension of the array B.  LDB >= max(1,N).
!
!    [out]   INFO  
!
!              INFO is INTEGER
!              = 0:  successful exit
!              < 0:  if INFO = -i, the i-th argument had an illegal value
!              > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                    has been completed, but the factor U is exactly
!                    singular, so the solution could not be computed.
      integer(kind=4) N
      integer(kind=4) IPIV(N)
      integer i, j
      integer INFO
      real(qp) C(N,N), B(N), U(N)
      real(dp) TEMP(N,N), U_TEMP(N)

      do i=1,N
         do j=1,N
            TEMP(i,j) = real( C(i,j), kind=qp)
         enddo
      enddo
      do i=1,N
         U_TEMP(i) = real( B(i), kind=qp)
         U(i) =  B(i)
      enddo
!      write(*,'(A)') '------------------------' 
!      do i=1,N
!      do j=1,N
!         write(*,'(e11.4, e11.4)') C(i,j), TEMP(i,j) 
!      enddo
!      enddo
      INFO = 0
      ! dgesv_(&eight, &one, A_cm, &eight, ipiv, b, &eight, &info);
      !    call matrix_out(u,8,1,002)
      write(*,'(A)') '------------------------' 
      do j=1,N
         write(*,'(e11.4)', advance='no') U_TEMP(j)
      enddo
      write(*,*)
      call dgesv(N, 1, TEMP, N, IPIV, U_TEMP, N, INFO )
      do j=1,N
         write(*,'(e11.4)', advance='no') U(j)
      enddo
      write(*,*)
      call qgesl(C, N, N, IPIV, u,0)

      if (INFO .ne. 0) then
         write(*,'(A,I2.1)') 'This is the value of INFO',INFO
         write(*,*) " Problem solving the system"
         write(*,*) "Exit!"
         stop
      endif
      do i=1,N
         U(i) = real( U_TEMP(i), kind=qp)
      enddo
      end subroutine solve_lgs

      subroutine matunitinit(A,n,m)
      use constants
      implicit none
      
      integer n, m, i, j
      real(qp) A(n,m)

      do i=1,n
         do j=1,m
            if (i .eq. j) then
               A(i,j) = 1.0
            else
               A(i,j) = 0.0
            endif
         enddo
      enddo

      end subroutine matunitinit
