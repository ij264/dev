   integer,       parameter :: i4b = selected_int_kind(9)
   ! Float and integers definition
   integer,       parameter :: dp  = selected_real_kind(8) 
   !> below was the definition from linpack routines
   integer,       parameter :: qp  = selected_real_kind(33,4931)

   !> math/physical constants
   real(qp),      parameter :: grv_uni_const  = 6.67e-11
   real(dp),      parameter :: pi=3.141592653589793238462643383279502884197
   real(dp),      parameter :: sqpi = 1.0/sqrt(pi)
   real(dp),      parameter :: Rearth=6370000.0
   real(qp),      parameter :: Emass=5.974e24
   !> factor for conversion of anomalies (Spherical Coefficients)
   real(dp),      parameter :: sqt=sqrt(4.0*pi)
   !> conversion of rad2grad or Vice Versa
   real(dp),      parameter :: deg2rad = pi/180
   real(dp),      parameter :: rad2deg = 180/pi

   !> ------- SIZE PARAMETERS -----------
   !> maximum spherical harmonic degree in the density file 
   integer(i4b),  parameter :: l_max = 50
   !> minimum spherical harmonic degree for computation 
   integer(i4b),  parameter :: l_min = 1
   !> maximum spherical harmonic degree for computation 
   integer(i4b),  parameter :: l_max_calc = 30
   !> number of radial layers (with thickness)
   !> with nr+1 number of radial grid points
   integer(i4b),  parameter :: nr = 256 
   !> half-width of the gaussian function
   !> for smoothing the viscosity profile
   integer(i4b),  parameter :: hlf_ny = 100
   !> Sigma for making the gaussian function
   real(qp), parameter :: gaus_sigma = 5.0

   !> ------- TYPE OF CALCULATION  -------
   !> flag deciding if we are doing an inversion 
   !> inv_flg ==1,2 : on, inv_flg == 0 just a fwd calculation
   integer(i4b),  parameter :: inv_flg = 0 

   !> ------- PARAMETERS FOR A POWELL -----
   !> ------- TYPE OF INVERSION  ----------
   !> ----------- INV_FLG==2 --------------
   !> Number of parameters to the inversion for
   !> for example, upper and lower mantle viscosity inversion
   integer(i4b),  parameter :: np=3
   real(qp),      dimension(np)  :: prms_begin
   real(qp),      dimension(np+1):: prms_bounds_dpth 
   real(qp),      dimension(np+1):: prms_bounds_rshl 
   integer(i4b),  dimension(np+1):: prms_bounds_indx 
   !> after how many iterations should parameters be written out
   integer(i4b),  parameter :: inv_out_skp=0
   !> here are the boundaries for the parameters

   !> ----------- GRID PARAMETERS -----------
   integer(i4b),  parameter :: grid = 179
   integer(i4b),  parameter :: grid2= 359
   real(dp),      parameter :: lat1 = -89
   real(dp),      parameter :: lat2 =  89 
   real(dp),      parameter :: lon1 =-179
   real(dp),      parameter :: lon2 = 179
   real(dp),      parameter :: Vspac=(lat2 - lat1)/(grid-1)
   real(dp),      parameter :: Hspac=(lon2 - lon1)/(grid2-1)

   !> ------------ OUTPUT PARAMETERS -------------------
   integer(i4b),  parameter:: kernel_output_flg= 1
   integer(i4b),  parameter:: visual3d = 0 
   !> if Correlation with RefGeo should be calculated
   integer(i4b),  parameter:: CorrFlag = 1
   !> ------ ACTUAL CALCULATION PARAMETERS -------------
   !> Which type kappaindex==1 compressible else: incompressible
   integer(i4b),  parameter:: kappaindex=1
   !> Starting Depth for calculating the gravity based on (050/-10)
   real(qp),      parameter:: top_depth =  050
   real(qp),      parameter:: bot_depth =  -10
   !> boundary condition on top 
   integer(i4b),  parameter:: bdry=1
   !> boundary condition on bot 
   integer(i4b),  parameter:: bdry_bot=1
   !> background density
   real(qp),      parameter:: bgrho = 4.0e3
   !> background Viscosity
   real(qp),      parameter:: visc0 = 1.   
   !> the changes of the Dynamic Topography at CMB surface
   !> has to do with rho_c
   real(qp),      parameter:: rho_c =  9.90e3
   !> the Dynamic topography on the surface 
   real(qp),      parameter:: rho_w = 1.03e3

   !------ FOR A TWO LAYER MODEL WITH MANY CALCULATIONS -
   !> minimum viscosity in the asthenosphere [Pa.s]
   real(qp),      parameter :: visc_lith = 1e23;
   !> thickness of the lithosphere [m]
   real(qp),      parameter :: thck_lith = 70e3;
   !> minimum viscosity in the asthenosphere [Pa.s] 
   real(qp),      parameter :: min_visc_asth = 1e18;
   !> lower mantle viscosity [Pa.s]
   real(qp),      parameter :: visc_mantle = 1e22;
   !> minimum thickness to try for the asthenosphere [m]
   real(qp),      parameter :: min_thck_asth = 11e3
   !> maximum thickness to try for the asthenosphere [m]
   real(qp),      parameter :: max_thck_asth = 1.2e6;
   !> the number of parameters to try between the
   !>    maximum and minimum viscosity
   integer(i4b),  parameter :: int_visc = 50 
   !> the number of parameters to try between the
   !>    maximum and minimum thickness for the lithosphere
   integer(i4b),  parameter :: int_thck = 70 

