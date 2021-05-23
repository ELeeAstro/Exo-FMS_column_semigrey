!!!
! Elspeth KH Lee - May 2021
! A simple program that emulates a column inside the Exo-FMS GCM
! This is useful for testing different RT solutions
! This version is for semi-grey radiation only
!
! Input parameters are set via the namelist (FMS_RC.nml) - hopefully paramaters are
! somewhat self explanatory. See the github readme for more information.
! Not guarenteed to be bug free
! NOTE: Indexing starts at 1, where 1 is the Top Of Atmosphere level or layer
!!!

program Exo_FMS_RC
  use, intrinsic :: iso_fortran_env
  use ts_isothermal_mod, only : ts_isothermal
  use ts_Toon_mod, only : ts_Toon
  use ts_Heng_mod, only : ts_Heng
  use ts_short_char_mod, only : ts_short_char
  use ts_Mendonca_mod, only : ts_Mendonca
  use k_Rosseland_mod, only : k_Ross_TK19
  use IC_mod, only : IC_profile
  use dry_conv_adj_mod, only : Ray_dry_adj
  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp

  integer :: n, i, k, u, j, inan
  integer :: nstep, nlay, nlev
  real(dp) :: t_step, t_tot
  real(dp) :: mu_z, Tirr, Tint, F0, Fint, pref, pu
  real(dp) :: tau_Vref, tau_IRref
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe
  real(dp), allocatable, dimension(:) :: k_Vl, k_IRl
  real(dp), allocatable, dimension(:) :: tau_Ve, tau_IRe, tau_IRl
  real(dp), allocatable, dimension(:) :: dT_rad, dT_conv, net_F

  real(dp) :: cp_air, grav, k_IR, k_V, kappa_air, Rd_air

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  real(dp) :: fl, AB

  integer :: ua, ub
  character(len=50) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a, b

  real(dp) :: start, finish

  character(len=50) :: ts_scheme, opac_scheme, adj_scheme

  integer :: u_nml

  namelist /FMS_RC_nml/ ts_scheme, opac_scheme, adj_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, AB, fl, &
          & iIC, corr

  !! Read input variables from namelist
  open(newunit=u_nml, file='FMS_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_RC_nml)
  close(u_nml)

  !! Number of layer edges (levels)
  nlev = nlay + 1

  !! Read in hybrid sigma grid values
  open(newunit=ua,file=trim(a_sh), action='read', status='old')
  open(newunit=ub,file=trim(b_sh), action='read', status='old')
  allocate(a(nlev),b(nlev))
  do k = 1, nlev
    read(ua,*) a(k)
    read(ub,*) b(k)
  end do
  close(ua); close(ub)

  !! Contruct pressure array [pa] at the levels using the hybrid sigma formula
  ! Reference surface pressure [pa] is pref
  allocate(pe(nlev))
  do k = 1, nlev
    pe(k) = a(k) + b(k)*pref
  end do
  pu = pe(1)

  !! Pressure at the layers
  allocate(pl(nlay),dpe(nlay))
  do k = 1, nlay
    dpe(k) = pe(k+1) - pe(k)
    pl(k) = dpe(k) / log(pe(k+1)/pe(k))
  end do

  !! Allocate other arrays we need
  allocate(Tl(nlay), dT_rad(nlay), dT_conv(nlay), net_F(nlev))
  allocate(tau_Ve(nlev),tau_IRe(nlev), k_Vl(nlay), k_IRl(nlay))

  if (ts_scheme == 'Heng') then
    allocate(tau_IRl(nlay))
  end if

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp

  F0 = sb * Tirr**4 ! Substellar point irradiation flux
  Fint = sb * Tint**4 ! Internal flux

  print*, 'Tint ', 'Tirr ', 'pref ', 'pu ', 'mu_z ', 'grav '
  print*, Tint, Tirr, pref/1e5_dp, pu/1e5_dp, mu_z, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_Vl(:) = k_V
  k_IRl(:) = k_IR

  select case(opac_scheme)
  case('Constant')
    tau_Vref = (k_V * pref) / grav
    tau_IRref = (k_IR * pref) / grav
    fl = 1.0_dp
  case('Heng')
    tau_Vref = (k_V * pref) / grav
    tau_IRref = (k_IR * pref) / grav / fl
  case default
  end select

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_Vl,k_IRl,Tint,mu_z,Tirr,grav,fl,Tl,prc)

  !! Print initial T-p profile
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i)
  end do
  print*, '--------------'


  !! Write out initial conditions
  open(newunit=u,file='FMS_RC_ic.out',action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_conv(:) = 0.0_dp

    select case(opac_scheme)
    case('Constant')
      tau_Ve(:) = tau_Vref * pe(:)/pref
      tau_IRe(:) = tau_IRref * pe(:)/pref
    case('Heng')
      tau_Ve(:) = tau_Vref * pe(:)/pref
      tau_IRe(:) = fl*tau_IRref*(pe(:)/pref)  + (1.0_dp - fl)*tau_IRref*(pe(:)/pref)**2
    case('TK19')
      ! Calculate optical depth structure for TK19 scheme
      tau_Ve(1) = 0.0_dp
      tau_IRe(1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_TK19(pl(k), k_Vl(k), k_IRl(k))
        tau_Ve(k+1) = tau_Ve(k) + (k_Vl(k) * dpe(k)) / grav
        tau_IRe(k+1) = tau_IRe(k) + (k_IRl(k) * dpe(k)) / grav
      end do
    case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select


    !! Two stream radiative transfer step
    select case(ts_scheme)
    case('Isothermal')
      ! Isothermal layers approximation
      call ts_isothermal(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, net_F)
    case('Toon')
      ! Toon method without scattering
      call ts_Toon(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, net_F)
    case('Shortchar')
      ! Short characteristics method without scattering
      call ts_short_char(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, net_F)
    case('Heng')
      ! Heng flux method without scattering
      tau_IRl(:) = fl*tau_IRref*(pl(:)/pref)  + (1.0_dp - fl)*tau_IRref*(pl(:)/pref)**2  ! Optical depth at layer midpoints
      call ts_Heng(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, tau_IRl, mu_z, F0, Tint, AB, net_F)
    case('Mendonca')
      !! In development !!
      ! Mendonca method without scattering
      call ts_Mendonca(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, net_F)
    case('None')
    case default
      print*, 'Invalid ts_scheme: ', trim(ts_scheme)
      stop
    end select

    !! Calculate the temperature tendency due to radiation
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    select case(adj_scheme)
    case('Ray_dry')
      ! Dry convective adjustment following Ray Pierrehumbert's python script
      call Ray_dry_adj(nlay, nlev, t_step, kappa_air, Tl, pl, pe, dT_conv)
    case('None')
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Check for NaN's in the temperature and exit the simulation if detected
    do i = 1, nlay
      if (ieee_is_nan(Tl(i)) .eqv. .True.) then
        do j = 1, nlay
          print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
        end do
        print*, nlev, net_F(nlev)
        inan = 1
        exit
      end if
    end do
    if (inan == 1) then
      exit
    end if

    !! Increase the total time simulated
    t_tot = t_tot + t_step

  end do

  !! cpu timer end
  call cpu_time(finish)

  !! Output the results
  print*, 'sec: ', 'hours: ', 'days: '
  print*, t_tot, t_tot/60.0_dp/60.0_dp,t_tot/60.0_dp/60.0_dp/24.0_dp

  print*, 'For profile properties: '
  print*, Tint, Tirr, pref, mu_z

  print*, 'Outputting results: '
  open(newunit=u,file='FMS_RC_pp.out', action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i), dT_rad(i), dT_conv(i), 0.5_dp*(tau_Ve(i+1)+tau_Ve(i)), 0.5_dp*(tau_IRe(i+1)+tau_IRe(i))
  end do
  close(u)

  print*, nstep, 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start


end program Exo_FMS_RC