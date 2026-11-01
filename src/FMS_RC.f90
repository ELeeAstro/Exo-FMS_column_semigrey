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
  ! sw_SH2 and sw_SH4 are disabled until their matrix assembly and flux
  ! reconstruction are repaired and validated.
  use sw_disort_ts_mod, only : sw_disort_ts
  use sw_disort_mod, only : sw_disort
  !use sw_Feautrier_mod, only : sw_Feautrier

  use lw_AA_E_mod, only : lw_AA_E
  use lw_AA_L_mod, only : lw_AA_L
  use lw_sc_linear_mod, only : lw_sc_linear
  !use lw_sc_scat_mod, only : lw_sc_scat   ! lw_sc_scat disabled (in development)
  use lw_sc_scat_ric_mod, only : lw_sc_scat_ric
  !use lw_sc_scat_iter_mod, only : lw_sc_scat_iter   ! lw_sc_scat_iter disabled (in development)
  use lw_VIM_E_mod, only : lw_VIM_E
  use lw_VIM_L_mod, only : lw_VIM_L
  use lw_Toon_mod, only : lw_Toon
  use lw_SH2_mod, only : lw_SH2
  use lw_SH4_mod, only : lw_SH4
  use lw_disort_ts_mod, only : lw_disort_ts
  use lw_disort_mod, only : lw_disort
  use lw_AD_mod, only : lw_AD
  !use lw_Feautrier_mod, only : lw_Feautrier   ! lw_Feautrier disabled (in development)
  use lw_DFE_mod, only : lw_DFE

  use k_Rosseland_mod, only : k_Ross_TK19, k_Ross_Freedman, k_Ross_Valencia

  use IC_mod, only : IC_profile

  use Manabe_dry_adjust_mod, only : Manabe_dry_adjust
  use Pierrehumbert_dry_adjust_mod, only : Pierrehumbert_dry_adjust
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

  real(dp), allocatable, dimension(:) :: Kzz, Kzz_e, Tl_rad, Rd_bar, cp_bar, kappa_bar

  real(dp) :: cp_air, grav, k_IR, k_V, kappa_air, Rd_air
  real(dp) :: sw_ac, sw_gc, lw_ac, lw_gc

  logical :: zcorr
  integer :: zcorr_meth
  real(dp) :: radius
  real(dp), allocatable, dimension(:) :: mu_z_eff, alt, alp

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  real(dp) :: fl, AB, olr, asr

  real(dp) :: Tdeep_init, Tstrat_init, ns, nl

  integer :: ua, ub
  character(len=50) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a, b

  real(dp) :: start, finish

  character(len=50) :: sw_scheme, lw_scheme, opac_scheme, adj_scheme
  character(len=20) :: lw_sc_iter_method
  integer :: lw_sc_iter_maxiter, lw_sc_iter_restart, lw_sc_iter_iterations
  real(dp) :: lw_sc_iter_tol, lw_sc_iter_residual

  integer :: u_nml

  namelist /FMS_RC_nml/ sw_scheme, lw_scheme, opac_scheme, adj_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, AB, fl, met, &
          & iIC, corr, sw_ac, sw_gc, lw_ac, lw_gc, Tdeep_init, Tstrat_init, ns, nl, &
          & zcorr, zcorr_meth, radius, lw_sc_iter_method, lw_sc_iter_tol, &
          & lw_sc_iter_maxiter, lw_sc_iter_restart

  !! Read input variables from namelist
  lw_sc_iter_method = 'ali'
  lw_sc_iter_tol = 1.0e-10_dp
  lw_sc_iter_maxiter = 10000
  lw_sc_iter_restart = 30
  lw_sc_iter_iterations = 0
  lw_sc_iter_residual = 0.0_dp
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
  allocate(Kzz(nlay), Kzz_e(nlev), Tl_rad(nlay), Rd_bar(nlay), cp_bar(nlay), kappa_bar(nlay))

  !! Allocate cloud properties - assume constant for all layers in this test code
  !! Adjust these parts for your own situations
  allocate(sw_a(nlay), sw_g(nlay), lw_a(nlay), lw_g(nlay))
  sw_a(:) = sw_ac
  sw_g(:) = sw_gc
  lw_a(:) = lw_ac
  lw_g(:) = lw_gc

  !! Deterministic initial values for all time-dependent fluxes and diagnostics.
  !! These also define valid output when a radiation or adjustment scheme is
  !! explicitly disabled, or when nstep is zero.
  dT_rad(:) = 0.0_dp
  dT_conv(:) = 0.0_dp
  net_F(:) = 0.0_dp
  tau_Ve(:) = 0.0_dp
  tau_IRe(:) = 0.0_dp
  sw_up(:) = 0.0_dp
  sw_down(:) = 0.0_dp
  sw_net(:) = 0.0_dp
  lw_up(:) = 0.0_dp
  lw_down(:) = 0.0_dp
  lw_net(:) = 0.0_dp
  Kzz(:) = 0.0_dp
  alt(:) = 0.0_dp
  mu_z_eff(:) = 0.0_dp
  alp(:) = 0.0_dp
  asr = 0.0_dp
  olr = 0.0_dp

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp

  Rd_bar(:) = Rd_air
  cp_bar(:) = cp_air
  kappa_bar(:) = kappa_air

  !! The prescribed Bond albedo removes a fixed fraction of the incident
  !! stellar flux before it enters the column. This is the usual semi-grey
  !! convention (e.g. Guillot 2010) and is intended for absorption-only
  !! shortwave schemes. A scattering sw_scheme computes the reflected flux
  !! explicitly, so leave AB = 0 there to avoid counting reflection twice.
  if ((AB < 0.0_dp) .or. (AB >= 1.0_dp)) then
    print*, 'Invalid AB (Bond albedo), must be in [0,1): ', AB
    stop
  end if

  Finc = (1.0_dp - AB) * sb * Tirr**4 ! Substellar point irradiation flux
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
    & Tdeep_init, Tstrat_init, kappa_air)

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
  print*, 'Using: ', trim(sw_scheme), ' + ' , trim(lw_scheme)
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0
  n = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_rad(:) = 0.0_dp
    dT_conv(:) = 0.0_dp
    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp
    sw_net(:) = 0.0_dp
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    lw_net(:) = 0.0_dp
    asr = 0.0_dp
    olr = 0.0_dp

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
      call sw_direct(nlev, tau_Ve, mu_z_eff, Finc, sw_up, sw_down, sw_net, asr)
    case('sw_adding')
      ! Approximate adding method with approximate scattering
      call sw_adding(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_SDA')
      ! Spherical harmonic doubling (SDA) adding four stream method with approximate scattering
      call sw_SDA(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_Toon')
      ! Toon89 shortwave method with multiple scattering
      call sw_Toon(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_disort_ts')
      ! Two stream disort method with multiple scattering
      call sw_disort_ts(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_disort')
      ! n-stream disort method with multiple scattering - benchmark method
      call sw_disort(nlay, nlev, tau_Ve, mu_z_eff, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_SH2','sw_SH4')
      print*, 'Unsupported sw_scheme (disabled due known correctness defects): ', trim(sw_scheme)
      error stop
    case('sw_Feautier','sw_Feautrier','sw_AD')
      print*, 'Unsupported sw_scheme (not implemented): ', trim(sw_scheme)
      error stop
    case('none')
      sw_up(:) = 0.0_dp
      sw_down(:) = 0.0_dp
      sw_net(:) = 0.0_dp
      asr = 0.0_dp
    case default
      print*, 'Invalid sw_scheme: ', trim(sw_scheme)
      stop
    end select

    !! A gas giant has no reflecting solid surface: stellar radiation that
    !! reaches the lower boundary leaves the modeled atmosphere downward.
    if (trim(sw_scheme) /= 'none') then
      if (abs(sw_up(nlev)) > 1.0e-10_dp*max(1.0_dp,abs(Finc))) then
        print*, 'Invalid shortwave lower-boundary upward flux: ', sw_up(nlev)
        error stop
      end if
    end if

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
      call lw_sc_linear(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)     
    case('lw_sc_scat')
      ! Short characteristics (sc) discrete ordinates with multiple scattering, direct block-tridiagonal solve
      ! - disabled (in development)
      !call lw_sc_scat(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_sc_scat_ric')
      ! Same SC equations solved by Riccati invariant imbedding
      call lw_sc_scat_ric(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_sc_scat_iter')
      ! Iterative SC: ALI, Ng-accelerated ALI, symmetric GS, or preconditioned GMRES
      ! - disabled (in development)
      !call lw_sc_scat_iter(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, &
      !  & lw_sc_iter_method, lw_sc_iter_tol, lw_sc_iter_maxiter, lw_sc_iter_restart, &
      !  & lw_up, lw_down, lw_net, olr, lw_sc_iter_iterations, lw_sc_iter_residual)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_VIM_E')
      ! Variational Iteration Method (VIM) with approximate scattering
      call lw_VIM_E(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_VIM_L')
      ! Variational Iteration Method (VIM), linear in tau Planck function
      call lw_VIM_L(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Toon')
      ! Toon89 longwave method with multiple scattering
      call lw_Toon(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_SH2')
      ! Two-stream spherical harmonic method (PICASO) with multiple scattering
      call lw_SH2(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_SH4')
      ! Four-stream spherical harmonic method (PICASO) with multiple scattering
      call lw_SH4(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort_ts')
      ! Two stream disort method with multiple scattering
      call lw_disort_ts(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort')
      ! n-stream disort method with multiple scattering - benchmark method
      call lw_disort(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_AD')
      ! Adding-doubling method with multiple scattering
      call lw_AD(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Feautrier')
      ! Feautrier method with ALI and Ng acceleration (petitRADTRANS), isotropic scattering - disabled (in development)
      !call lw_Feautrier(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_DFE')
      ! Absorption-only Discontinuous Finite Element (DFE) method
      if (any(abs(lw_a(:)) > 1.0e-12_dp)) then
        print*, 'lw_DFE is absorption-only; set lw_ac = 0.'
        error stop
      end if
      call lw_DFE(nlay, nlev, Tl, pl, pe, tau_IRe, Tint, lw_up, lw_down, lw_net, olr)
    case('none')
      lw_up(:) = 0.0_dp
      lw_down(:) = 0.0_dp
      lw_net(:) = 0.0_dp
      olr = 0.0_dp
    case default
      print*, 'Invalid lw_scheme: ', trim(lw_scheme)
      stop
    end select

    !! Every gas-giant longwave solver must expose the prescribed internal
    !! net flux at the lower boundary.
    if (trim(lw_scheme) /= 'none') then
      ! The legacy DISORT solver closes this condition to about 1e-8
      ! relative precision because its band-integrated Planck function is
      ! evaluated numerically.
      if (abs(lw_net(nlev) - Fint) > 1.0e-7_dp*max(1.0_dp,abs(Fint))) then
        print*, 'Invalid longwave lower-boundary flux [actual, expected]: '
        print*, lw_net(nlev), Fint
        error stop
      end if
    end if

    !! Calculate the temperature tendency due to radiation
    net_F(:) = lw_net(:) + sw_net(:) ! Net flux into/out of level
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    !! The schemes act on the column after this step's radiative heating,
    !! so the column is convectively stable at the end of the step
    Tl_rad(:) = Tl(:) + t_step * dT_rad(:)
    select case(adj_scheme)
    case('Manabe')
      ! Dry convective adjustment after Manabe, Smagorinsky & Strickler (1965), no iteration
      call Manabe_dry_adjust(nlay, Tl_rad, pl, pe, kappa_air, t_step, dT_conv)
      Kzz(:) = 1e1_dp
    case('Pierrehumbert','Ray_dry')
      ! Pairwise dry convective adjustment following Ray Pierrehumbert's python script
      call Pierrehumbert_dry_adjust(nlay, Tl_rad, pl, pe, kappa_air, t_step, dT_conv)
      Kzz(:) = 1e1_dp
    case('MLT')
      ! Use mixing length theory (MLT) to time dependently adjust the adiabat and estimate Kzz
      call MLT(nlay, Tl_rad, pl, pe, Rd_bar, cp_bar, grav, t_step, dT_conv, Kzz=Kzz_e)
      ! MLT gives Kzz at the levels, output the layer average
      Kzz(:) = 0.5_dp * (Kzz_e(1:nlay) + Kzz_e(2:nlev))
    case('none')
      dT_conv(:) = 0.0_dp
      Kzz(:) = 1e1_dp
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
        print*, nlev, net_F(nlev), olr
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
  print*, Fint, Tint

  if (trim(lw_scheme) == 'lw_sc_scat_iter') then
    print*, 'Iterative SC method, iterations, relative residual:'
    print*, trim(lw_sc_iter_method), lw_sc_iter_iterations, lw_sc_iter_residual
  end if

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
