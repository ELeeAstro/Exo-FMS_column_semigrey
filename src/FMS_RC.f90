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

  use sw_direct_mod, only : sw_direct
  use sw_adding_mod, only : sw_adding
  use sw_SDA_mod, only : sw_SDA
  use sw_Toon_mod, only : sw_Toon
  use sw_SH2_mod, only : sw_SH2
  use sw_SH4_mod, only : sw_SH4
  use sw_disort_ts_mod, only : sw_disort_ts
  !use sw_disort_mod, only : sw_disort
  !use sw_Feautrier_mod, only : sw_Feautrier

  use lw_AA_E_mod, only : lw_AA_E
  use lw_AA_L_mod, only : lw_AA_L
  use lw_sc_linear_mod, only : lw_sc_linear
  use lw_VIM_mod, only : lw_VIM
  use lw_Toon_mod, only : lw_Toon
  use lw_disort_ts_mod, only : lw_disort_ts
  !use lw_disort_mod, only : lw_disort
  use lw_Feautrier_mod, only : lw_Feautrier
  !use lw_AD_mod, only : lw_AD
  use lw_DFE_mod, only : lw_DFE

  use k_Rosseland_mod, only : k_Ross_TK19, k_Ross_Freedman, k_Ross_Valencia

  use IC_mod, only : IC_profile

  use dry_conv_adj_mod, only : Ray_dry_adj
  use MLT_mod, only : MLT

  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp

  integer :: n, i, k, u, uu, j, inan
  integer :: nstep, nlay, nlev
  real(dp) :: t_step, t_tot
  real(dp) :: mu_z, Tirr, Tint, Finc, Fint, pref, pu, met
  real(dp) :: tau_Vref, tau_IRref
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe
  real(dp), allocatable, dimension(:) :: k_Vl, k_IRl
  real(dp), allocatable, dimension(:) :: tau_Ve, tau_IRe
  real(dp), allocatable, dimension(:) :: dT_rad, dT_conv, net_F

  real(dp), allocatable, dimension(:) :: sw_a, sw_g, lw_a, lw_g

  real(dp), allocatable, dimension(:) :: sw_up, sw_down, sw_net, lw_up, lw_down, lw_net

  real(dp), allocatable, dimension(:) :: Kzz, Rd_bar, cp_bar, kappa_bar

  real(dp) :: cp_air, grav, k_IR, k_V, kappa_air, Rd_air
  real(dp) :: sw_ac, sw_gc, lw_ac, lw_gc

  logical :: zcorr
  integer :: zcorr_meth
  real(dp) :: radius
  real(dp), allocatable, dimension(:) :: mu_z_eff, alt, alp

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  real(dp) :: fl, AB, olr, asr, Ts, dTs, net_Fs

  logical :: surf
  real(dp) :: Ts_init, sw_a_surf, cp_surf, Tstrat_init, lw_a_surf, ns, nl

  integer :: ua, ub
  character(len=50) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a, b

  real(dp) :: start, finish

  character(len=50) :: sw_scheme, lw_scheme, opac_scheme, adj_scheme

  integer :: u_nml

  namelist /FMS_RC_nml/ sw_scheme, lw_scheme, opac_scheme, adj_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, AB, fl, met, &
          & iIC, corr, sw_ac, sw_gc, lw_ac, lw_gc, sw_a_surf, surf, Ts_init, cp_surf, &
          & Tstrat_init, lw_a_surf, ns, nl, zcorr, zcorr_meth, radius

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
  allocate(alt(nlev), mu_z_eff(nlev), alp(nlev))
  allocate(sw_up(nlev), sw_down(nlev), sw_net(nlev), lw_up(nlev), lw_down(nlev), lw_net(nlev))
  allocate(Kzz(nlay), Rd_bar(nlay), cp_bar(nlay), kappa_bar(nlay))

  !! Allocate cloud properties - assume constant for all layers in this test code
  !! Adjust these parts for your own situations
  allocate(sw_a(nlay), sw_g(nlay), lw_a(nlay), lw_g(nlay))
  sw_a(:) = sw_ac
  sw_g(:) = sw_gc
  lw_a(:) = lw_ac
  lw_g(:) = lw_gc

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp

  Rd_bar(:) = Rd_air
  cp_bar(:) = cp_air
  kappa_bar(:) = kappa_air

  Finc = sb * Tirr**4 ! Substellar point irradiation flux
  Fint = sb * Tint**4 ! Internal flux

  print*, 'Tint ', 'Tirr ', 'pref ', 'pu ', 'mu_z ', 'grav '
  print*, Tint, Tirr, pref/1e5_dp, pu/1e5_dp, mu_z, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_Vl(:) = k_V
  k_IRl(:) = k_IR

  select case(opac_scheme)
  case('constant')
    tau_Vref = (k_V * pref) / grav
    tau_IRref = (k_IR * pref) / grav
    fl = 1.0_dp
  case('Heng')
    tau_Vref = (k_V * pref) / grav
    tau_IRref = (k_IR * pref) / grav / fl
  case default
    tau_Vref = (k_V * pref) / grav
    tau_IRref = (k_IR * pref) / grav
    fl = 1.0_dp
  end select

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_Vl,k_IRl,Tint,mu_z,Tirr,grav,fl,Tl,prc, &
    & Ts_init, Tstrat_init, kappa_air)

  !! Print initial T-p profile
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i)
  end do
  print*, '--------------'

  Ts = Ts_init

  !! Write out initial conditions
  open(newunit=u,file='FMS_RC_ic.out',action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Using: ', trim(sw_scheme), ' + ' , trim(lw_scheme)
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_conv(:) = 0.0_dp

    select case(opac_scheme)
    case('constant')
      tau_Ve(:) = tau_Vref * pe(:)/pref
      tau_IRe(:) = tau_IRref * pe(:)/pref
    case('Heng')
      tau_Ve(:) = tau_Vref * (pe(:)/pref)**ns
      tau_IRe(:) = fl*tau_IRref*(pe(:)/pref)  + (1.0_dp - fl)*tau_IRref*(pe(:)/pref)**nl
    case('TK19')
      ! Calculate optical depth structure for TK19 scheme
      tau_Ve(1) = 0.0_dp
      tau_IRe(1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_TK19(pl(k), k_Vl(k), k_IRl(k))
        tau_Ve(k+1) = tau_Ve(k) + (k_Vl(k) * dpe(k)) / grav
        tau_IRe(k+1) = tau_IRe(k) + (k_IRl(k) * dpe(k)) / grav
      end do
    case('Freedman')
      ! Calculate optical depth structure for Freedman et al. (2014) Rosseland mean fitting function scheme
      tau_Ve(:) = tau_Vref * pe(:)/pref
      tau_IRe(1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_Freedman(Tl(k), pl(k), met, k_IRl(k))
        tau_IRe(k+1) = tau_IRe(k) + (k_IRl(k) * dpe(k)) / grav
      end do
    case('Valencia')
      ! Calculate optical depth structure for Valencia et al. (2013) Rosseland mean fitting function scheme
      tau_Ve(:) = tau_Vref * pe(:)/pref
      tau_IRe(1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_Valencia(Tl(k), pl(k), met, k_IRl(k))
        tau_IRe(k+1) = tau_IRe(k) + (k_IRl(k) * dpe(k)) / grav
      end do
    case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select


    !! Zenith angle geometric correction
    if ((zcorr .eqv. .True.) .and. (mu_z > 0.0_dp)) then
      ! First calculate the altitude at each level from the hypsometric equation
      ! Assume constant gravity for simplicity
      alt(nlev) = 0.0_dp
      do k = nlay, 1, -1
        alt(k) = alt(k+1) + (Rd_air*Tl(k))/grav * log(pe(k+1)/pe(k))
      end do

      select case(zcorr_meth)
      case(1)
        ! Basic geometric correction following Li & Shibata (2006) Eq. (2)
        mu_z_eff(:) = sqrt(1.0 - (radius/(radius + alt(:)))**2 * (1.0 - mu_z**2))
      case(2)
        ! Spherical layer correction following Li & Shibata (2006) Eq.(10)
        alp(nlev) = (alt(nlay) -  alt(nlev))/radius
        do k = nlay,1,-1
           alp(k) = (alt(k) -  alt(k+1))/(radius + alt(k))
        end do
        mu_z_eff(:) = (sqrt(1.0 - (radius/(radius + alt(:)))**2 * (1.0 - mu_z**2)) + &
        & sqrt((1.0 + alp(:))**2 - (radius/(radius + alt(:)))**2 * (1.0 - mu_z**2))) / &
        & (2.0 + alp(:))
      case default
        print*, 'Invalid zcorr_meth ', zcorr_meth
        stop
      end select
    else
      ! No correction, use single zenith angle
      mu_z_eff(:) = mu_z
    end if

    !! Shortwave radiative transfer step
    select case(sw_scheme)
    case('sw_direct')
      ! Direct beam only method with no scattering
      call sw_direct(nlev, tau_Ve, mu_z_eff, Finc, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_adding')
      ! Approximate adding method with approximate scattering
      call sw_adding(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SDA')
      ! Spherical harmonic doubling (SDA) adding four stream method with approximate scattering
      call sw_SDA(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_Toon')
      ! Toon89 shortwave method with multiple scattering
      call sw_Toon(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SH2')
      ! Spherical harmonic two stream method with multiple scattering
      call sw_SH2(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SH4')
      ! Spherical harmonic four stream method with multiple scattering
      call sw_SH4(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_disort_ts')
      ! Two stream disort method with multiple scattering
      call sw_disort_ts(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_disort')
      ! n-stream disort method with multiple scattering (benchmark method)
      !call sw_disort(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    case('sw_Feautier')
      ! Feautier method
      !call sw_Feautier(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)  
    case('sw_AD')
      ! Adding-doubling method
      !call sw_AD(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    case('none')
    case default
      print*, 'Invalid sw_scheme: ', trim(sw_scheme)
      stop
    end select

    !! Longwave radiative transfer step
    select case(lw_scheme)
    case('lw_AA_E')
      ! Absorption approximation exponential in tau (AA_E) with approximate scattering
      call lw_AA_E(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr) 
    case('lw_AA_L')
      ! Absorption approximation linear in tau (AA_L) with approximate scattering
      call lw_AA_L(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)     
    case('lw_sc_linear')
      ! Short characteristics (sc) with linear interpolants with no scattering
      call lw_sc_linear(nlay, nlev, Tl, pl, pe, tau_IRe, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_VIM')
      ! Variational Iteration Method (VIM) with approximate scattering
      call lw_VIM(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Toon')
      ! Toon89 longwave method with multiple scattering
      call lw_Toon(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, lw_a_surf, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort_ts')
      ! Two stream disort method with multiple scattering
      call lw_disort_ts(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort')
      ! n-stream disort method with multiple scattering (benchmark method)
      !call lw_disort(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr) 
    case('lw_Feautrier')
      ! Feautrier method
      call lw_Feautrier(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g ,Tint, lw_up, lw_down, lw_net, olr) 
    case('lw_AD')
      ! Adding-doubling method
      !call lw_AD(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, a_surf, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_DFE')
      ! Discontinuous Finite Element (DFE) method
      call lw_DFE(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g ,Tint, lw_up, lw_down, lw_net, olr)
    case('none')
    case default
      print*, 'Invalid lw_scheme: ', trim(lw_scheme)
      stop
    end select

    !! Calculate the temperature tendency due to radiation
    net_F(:) = lw_net(:) + sw_net(:) ! Net flux into/out of level
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    select case(adj_scheme)
    case('Ray_dry')
      ! Dry convective adjustment following Ray Pierrehumbert's python script
      call Ray_dry_adj(nlay, nlev, t_step, kappa_air, Tl, pl, pe, dT_conv)
      Kzz(:) = 1e1_dp
    case('MLT')
      ! Use mixing length theory (MLT) to time dependently adjust the adiabat and estimate Kzz
      call MLT(nlay, nlev, t_step, Tl, pl, pe, Rd_bar, cp_bar, kappa_bar, &
         & grav, dT_conv, Kzz)      
    case('none')
      Kzz(:) = 1e1_dp
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Forward march the temperature change in the surface
    if (surf .eqv. .True.) then
      dTs = net_Fs / (cp_surf)
      Ts = Ts + t_step * dTs
    else
      Ts = 0.0_dp
    end if

    !! Check for NaN's in the temperature and exit the simulation if detected
    do i = 1, nlay
      if (ieee_is_nan(Tl(i)) .eqv. .True.) then
        do j = 1, nlay
          print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
        end do
        print*, nlev, net_F(nlev), olr, net_Fs, Ts
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

  print*, 'OLR [W m-2], Teff:'
  print*, olr, (olr/sb)**(0.25_dp)

  print*, 'ASR [W m-2], Tinc:'
  print*, asr, (asr/sb)**(0.25_dp)

  print*, 'Internal T [W m-2], Tint'
  print*, sb * Tint**4, Tint

  print*, 'Surface T [K]: '
  print*, Ts

  print*, 'Outputting results: '
  open(newunit=uu,file='FMS_RC_pp.out', action='readwrite')
  do i = 1, nlay
    write(uu,*) i, pl(i), Tl(i), dT_rad(i), dT_conv(i), Kzz(i)
  end do
  close(uu)

  open(newunit=uu,file='FMS_RC_flx.out', action='readwrite')
  do i = 1, nlev
    write(uu,*) i, pe(i), sw_up(i), sw_down(i), sw_net(i), lw_up(i), lw_down(i), lw_net(i), & 
      & tau_Ve(i), tau_IRe(i)
  end do
  close(uu)

  print*, n, 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start

end program Exo_FMS_RC
