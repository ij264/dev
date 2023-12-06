      subroutine fldcalc(GeoK, GravK, VelK, SurK, CMBK)

      use constants
      use sph_tools, only : theta, phi, facto, leg 
      use ref_flds, only : georfCO, grvrfCO, dynrfCO, cmbrfCO
      use util, only : r, rho0, anom, sphcout, fldout
      use commpi

      implicit none
      
      common /path/ output, input
      common /ResCorr/ GeoCorr,GravCorr,DyntCorr,CMBCorr, surmax, geomax
      common /ResStat/ flatness, obs_surmax, obs_geomax 

      character*120 output, input

      real(dp) sqr, powone
      real(dp) geomax, surmax
      real(dp) sphharm
      ! Kernels Geo, Vel, Surf, CMB, GravK
      real(dp) GeoK  (l_min:l_max_calc, nr+1)
      real(dp) VelK  (l_min:l_max_calc, nr+1)
      real(dp) SurK  (l_min:l_max_calc, nr+1)
      real(dp) CMBK  (l_min:l_max_calc, nr+1)
      real(dp) GravK (l_min:l_max_calc, nr+1)
      ! SPH coeffs  geoid
      real(dp) geoidCO  (l_min:l_max_calc, -l_max_calc:+l_max_calc)
      ! SPH coeffs  gravity
      real(dp) gravCO   (l_min:l_max_calc, -l_max_calc:+l_max_calc) 
      ! SPH coeffs  surface velocity
      real(dp) SurTopCO (l_min:l_max_calc, -l_max_calc:+l_max_calc)
      ! SPH coeffs  surface dynamic topography
      real(dp) CMBTopCO (l_min:l_max_calc, -l_max_calc:+l_max_calc)
      ! SPH coeffs  cmb dynamic topography      
      real(dp) veloCO   (l_min:l_max_calc, -l_max_calc:+l_max_calc) 

      real(dp) geo_grd  (grid , grid2) 
      real(dp) grv_grd  (grid , grid2) 
      real(dp) srfdt_grd(grid , grid2) 
      real(dp) cmbdt_grd(grid , grid2) 
      real(dp) srfvl_grd(grid , grid2)
      real(dp) georf_grd(grid , grid2)
      real(dp) grvrf_grd(grid , grid2)
      real(dp) dynrf_grd(grid , grid2)
      real(dp) cmbrf_grd(grid , grid2)
      real(dp) GeoCorr, GravCorr, DyntCorr, flatness, CMBCorr
      real(dp) start, finish, obs_surmax, obs_geomax
      real(dp) laysize
      integer i, l, m, j, proc_caller

      laysize = real((r(nr+1)- r(1))/(nr),kind=dp)

   
      geoidCO (:,:)  = 0.0
      gravCO  (:,:)  = 0.0
      veloCO  (:,:)  = 0.0
      SurTopCO(:,:)  = 0.0
      CMBTopCO(:,:)  = 0.0

      geo_grd  (:,:) = 0.0
      grv_grd  (:,:) = 0.0
      srfvl_grd(:,:) = 0.0
      srfdt_grd(:,:) = 0.0
      cmbdt_grd(:,:) = 0.0

      georf_grd(:,:) = 0.0
      grvrf_grd(:,:) = 0.0
      dynrf_grd(:,:) = 0.0
      cmbrf_grd(:,:) = 0.0

      GeoCorr        = 0.0
      GravCorr       = 0.0
      DyntCorr       = 0.0
      CMBCorr        = 0.0
      surmax         = 0.0
      geomax         = 0.0

      proc_caller = 0
      do l=l_min, l_max_calc
         do m= -l, l
            if (mod(proc_caller,nproc) == myproc) then
               do j=nr+1,1,-1
                  geoidCO(l,m) =
     &               geoidCO(l,m) + anom(l,m,j)*GeoK(l,j)*laysize
                  gravCO (l,m) =
     &               gravCO(l,m)  + anom(l,m,j)*GravK(l,j)*laysize
                  veloCO(l,m)  =
     &               veloCO(l,m)  + anom(l,m,j)*VelK(l,j)*laysize
                  SurTopCO(l,m)=
     &               SurTopCO(l,m)+ anom(l,m,j)*SurK(l,j)*laysize
                  CMBTopCO(l,m)=
     &               CMBTopCO(l,m)+ anom(l,m,j)*CMBK(l,j)*laysize
               enddo
               do i=1,grid
               do j=1,grid2
               if(m==0) then
                  sphharm = sqr(2*l+1)*0.5*sqpi*leg(l,m,i)
               elseif (m<0) then
                  sphharm = sqr(2*l+1)*sqpi/sqr(2)*facto(l+m)/facto(l-m)
     &                           *leg(l,-m,i)*cos(m*phi(j))
               else 
                  sphharm = sqr(2*l+1)*sqpi/sqr(2)*facto(l-m)/facto(l+m)
     &                           *leg(l,m,i)*sin(m*phi(j))
               endif
               geo_grd  (i,j) = geo_grd  (i,j) + geoidCO (l,m)*sphharm
               grv_grd  (i,j) = grv_grd  (i,j) + gravCO  (l,m)*sphharm
               srfvl_grd(i,j) = srfvl_grd(i,j) + veloCO  (l,m)*sphharm
               srfdt_grd(i,j) = srfdt_grd(i,j) + SurTopCO(l,m)*sphharm
               cmbdt_grd(i,j) = cmbdt_grd(i,j) + CMBTopCO(l,m)*sphharm
               georf_grd(i,j) = georf_grd(i,j) + georfCO (l,m)*sphharm
               grvrf_grd(i,j) = grvrf_grd(i,j) + grvrfCO (l,m)*sphharm
               dynrf_grd(i,j) = dynrf_grd(i,j) + dynrfCO (l,m)*sphharm
               cmbrf_grd(i,j) = cmbrf_grd(i,j) + cmbrfCO (l,m)*sphharm
               enddo
               enddo
            endif
            proc_caller = proc_caller +1 
         enddo
      enddo

      call mpi_synch_coef(geoidCO)
      call mpi_synch_coef(gravCO )
      call mpi_synch_coef(veloCO )
      call mpi_synch_coef(SurTopCO)
      call mpi_synch_coef(CMBTopCO)

      call mpi_synch_fld(geo_grd)
      call mpi_synch_fld(grv_grd)
      call mpi_synch_fld(srfvl_grd)
      call mpi_synch_fld(srfdt_grd)
      call mpi_synch_fld(cmbdt_grd)
      call mpi_synch_fld(georf_grd)
      call mpi_synch_fld(grvrf_grd)
      call mpi_synch_fld(dynrf_grd)
      call mpi_synch_fld(cmbrf_grd)
      
      if (CorrFlag ==1) then
         GeoCorr  =   sum(geoidCO(:,:)*georfCO(:,:))/
     &    (sqrt(sum(geoidCO**2)*sum(georfCO**2)))
         GravCorr =   sum(gravCO(:,:)*grvrfCO(:,:))/
     &    (sqrt(sum(gravCO**2)*sum(grvrfCO**2)))
         DyntCorr = sum(dynrfCO(:,:)*SurTopCO(:,:))/
     &    (sqrt(sum(dynrfCO**2)*sum(SurTopCO**2)))
         CMBCorr = sum(cmbrfCO(:,:)*CMBTopCO(:,:))/
     &    (sqrt(sum(cmbrfCO**2)*sum(CMBTopCO**2)))
         surmax = maxval(abs(srfdt_grd))
         geomax = maxval(abs(geo_grd))
         obs_surmax = maxval(abs(dynrf_grd))
         obs_geomax = maxval(abs(georf_grd))
         !flatness  =
!     &   ! sqrt(sum(SurTopCO(2:3,:)**2)/sum(SurTopCO(8: 9,:)**2))
      endif

      if (myproc == root_proc) then
      if (inv_flg==0) then      
         ! Output of Spherical Harmonics
         call sphcout(geoidCO,  'SPH_geoid'            )
         call sphcout(grvrfCO,  'SPH_gravity'          )
         call sphcout(SurTopCO, 'SPH_surftopo'         )
         call sphcout(CMBTopCO, 'SPH_cmbtopo'          )
         call sphcout(veloCO,   'SPH_surfvel'          )
         call sphcout(georfCO,  'SPH_REF_geoid'        )
         call sphcout(grvrfCO,  'SPH_REF_gravity'      )
         call sphcout(dynrfCO,  'SPH_REF_dyntopography')

         ! Output of Grid
         call fldout(geo_grd  , 'GRD_geoid'            )
         call fldout(grv_grd  , 'GRD_gravity'          )
         call fldout(srfdt_grd, 'GRD_surftopo'         )
         call fldout(cmbdt_grd, 'GRD_cmbtopo'          )
         call fldout(srfvl_grd, 'GRD_surfvel'          )
         call fldout(georf_grd, 'GRD_REF_geoid'        )     
         call fldout(grvrf_grd, 'GRD_REF_gravity'      )     
         call fldout(dynrf_grd, 'GRD_REF_dyntopography')     
         call fldout(cmbrf_grd, 'GRD_REF_cmbtopography')     
      endif
      endif
      call wait4others
      end subroutine


* function sqri
!     Needed for some of the calculations
      function sqr(i) result(sqr_val)
      use constants
      implicit none
      integer i
      real(dp) sqr_val
      sqr_val = sqrt(real(i))
      end function sqr

*  function powone
!  To reduce the load of (-1)^m
      function powone(m) result(pow_one)
      use constants
      implicit none
      integer m
      real(dp) pow_one
      if (mod(m,2) .eq. 0) then
         pow_one = 1.0
      else
         pow_one = -1.0
      endif
      return
      end function powone
