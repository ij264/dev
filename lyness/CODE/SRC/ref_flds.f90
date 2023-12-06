module ref_flds
   use constants
   use util, only : get_file_unit
   implicit none

   public :: ref_flds_init
   
   real(dp), protected, dimension(l_min:l_max_calc, -l_max_calc:+l_max_calc), save :: georfCO 
   real(dp), protected, dimension(l_min:l_max_calc, -l_max_calc:+l_max_calc), save :: grvrfCO 
   real(dp), protected, dimension(l_min:l_max_calc, -l_max_calc:+l_max_calc), save :: dynrfCO 
   real(dp), protected, dimension(l_min:l_max_calc, -l_max_calc:+l_max_calc), save :: cmbrfCO 

contains

   ! read in the reference fields, the reference geo_grd is 
   ! the total geo_grd, including the the hydrostatic shape
   ! of the Earth. We convert the hydrostatic coefficients
   ! to the dynamic ones by using the hydrostatic shape
   ! coeffs of Nakiboglu 1982 

   subroutine ref_flds_init()
   use constants
   implicit none

   integer  :: l, m, logunit
   integer  :: ierr, iostatus
   character(len=50) :: line_tmp
   real(dp) :: cnm, snm

   logunit = get_file_unit(999)
   open(logunit, file='./INPUT/RefGeo/eigen5c2'&
          ,action='read', status='old', iostat=ierr)
   do l=0, l_max_calc
      do m=0,l
         if (m==0) then
            if (l>= l_min) then
            read(logunit,*) georfCO(l,m)
            else
            read(logunit,*)
            endif
            if (l==2) then
               georfCO(l,m) = georfCO(l,m) + 1071.2e-6/sqrt(2.0*2.0+1.0) ! Chambat et al. 2010
            else if (l==4) then
               georfCO(l,m) = georfCO(l,m) -    2.96e-6/sqrt(2.0*4.0+1.0) ! Chambat et al. 2010
            endif
!             if (l==2) then
!                georfCO(l,m) = georfCO(l,m) + 1072.618e-6/sqrt(2.0*2.0+1.0) ! Nakiboglu 1982
!             else if (l==4) then
!                georfCO(l,m) = georfCO(l,m) -    2.992e-6/sqrt(2.0*4.0+1.0) ! Nakiboglu 1982
!             endif
         else
            if (l>= l_min) then
            read(logunit,*) georfCO(l,-m), georfCO(l,+m)
            else
            read(logunit,*)
            endif
         endif
      enddo
   enddo
   close(logunit)
   logunit = -1


   grvrfCO(l_min:l_max_calc,-l_max_calc:+l_max_calc) = &
            reshape([([(georfCO(l,m)*(l-1), l=l_min,l_max_calc )],&
            m=-l_max_calc,l_max_calc)],(/l_max_calc-l_min+1,2*l_max_calc+1/))

   georfCO  = Rearth*georfCO*sqt
   grvrfCO = 1e5*((grv_uni_const*real(Emass,kind=dp))/Rearth**2)*grvrfCO*sqt

   logunit = get_file_unit(999)
   open(logunit,&
            file='./INPUT/RefGeo/Hoggard_et_al_2017_reformulated.gfc',&
            action='read',status='old')
   dynrfCO(:,:) = 0.0
   !do l=0,l_max_calc !50 for when l_max_calc > 50
   do l=0,l_max_calc
       do m=0,l
           if (m==0) then
               if (l>= l_min) then
               read(logunit,*) dynrfCO(l,m)
               else
               read(logunit,*)
               endif
           else
               if (l>= l_min) then
               read(logunit,*) dynrfCO(l,-m), dynrfCO(l,+m)
               else
               read(logunit,*)
               endif
           endif
       enddo
   enddo
   close(logunit)
   logunit = -1
   dynrfCO = dynrfCO*1e3 ! converting kilometer to meter

   cmbrfCO(:,:) = 0.0
   logunit = get_file_unit(999)
   open(logunit,&
            file='./INPUT/RefGeo/Koelemeijer_2017_cmb.sph',&
            action='read',status='old')
   read(logunit,'(A)', IOSTAT=iostatus)  line_tmp
   do while (iostatus>=0) 
      if (line_tmp(1:5)=='SPHCO') then
         read(line_tmp(6:),*) l, m ,cnm, snm
         if (l .ge. l_min .and. l .le. l_max_calc) then
         if (m==0) then
            cmbrfCO(l,0) = cnm
         else
            cmbrfCO(l,+m) = snm
            cmbrfCO(l,-m) = cnm
         endif
         endif
      endif
      read(logunit,'(A)', IOSTAT=iostatus) line_tmp
   enddo
   close(logunit)
   logunit = -1
   end subroutine ref_flds_init

end module ref_flds
