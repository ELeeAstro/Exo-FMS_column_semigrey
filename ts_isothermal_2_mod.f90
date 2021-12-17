!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Dec 2021 : adding method
! sw: Adding layer method with scattering
! lw: Two-stream method following the isothermal layer approximation
!      Pros: Very fast, better at high optical depths than 1st isothermal method
!      Cons: Worse at low optical depth, worse at lower internal temperatures
!!!

module ts_isothermal_2_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  public :: ts_isothermal_2
  private :: lw_grey_updown, sw_grey_updown_adding

contains

  subroutine ts_isothermal_2(surf, nlay, nlev, Ts,  Tl, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, sw_a_surf, lw_a_surf, net_F, olr, asr, net_Fs)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB, sw_a_surf, lw_a_surf, Ts
    real(dp), dimension(nlay), intent(in) :: Tl
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g


    !! Output variables
    real(dp), intent(out) :: olr, asr, net_Fs
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlay) :: bl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve(:), mu_z, sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    bl(:) = sb * Tl(:)**4  ! Integrated planck function flux at layers
    if (surf .eqv. .True.) then
      be_int = sb * Ts**4
    else
      be_int = sb * Tint**4 ! Integrated planck function intensity for internal temperature
    end if
    call lw_grey_updown(surf, nlay, nlev, bl, be_int, tau_IRe(:), lw_a_surf, lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Net surface flux (for surface temperature evolution)
    !! We have to define positive as downward (heating) and cooling (upward) in this case
    net_Fs = sw_down(nlev) + lw_down(nlev) - lw_up(nlev)

    !! Outut olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_isothermal_2

  subroutine lw_grey_updown(surf, nlay, nlev, bl, be_int, tau_IRe, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: be_int, lw_a_surf
    real(dp), dimension(nlev), intent(in) :: tau_IRe
    real(dp), dimension(nlay), intent(in) :: bl

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k
    real(dp), dimension(nlay) :: dtau

    !! Prepare loop
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
    end do

    !! First do the downward loop
    lw_down(1) = 0.0_dp
    do k = 1, nlay
       lw_down(k+1) = (2.0_dp * bl(k) * dtau(k))/(2.0_dp + dtau(k)) + &
       & lw_down(k) * (2.0_dp - dtau(k))/(2.0_dp + dtau(k))
    end do

    !! Perform upward loop
    if (surf .eqv. .True.) then
      ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
      lw_up(nlev) = lw_down(nlev)*lw_a_surf + be_int
    else
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here the lw_a_surf is assumed to be = 1 as per the definition
      ! here we use the same condition but use intensity units to be consistent
      lw_up(nlev) = lw_down(nlev) + be_int
    end if
    do k = nlay, 1, -1
      lw_up(k) = (2.0_dp * bl(k) * dtau(k))/(2.0_dp + dtau(k)) + &
      & lw_up(k+1) * (2.0_dp - dtau(k))/(2.0_dp + dtau(k))
    end do


  end subroutine lw_grey_updown

  subroutine sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables
    integer :: k
    real(dp) :: lamtau, e_lamtau, lim, arg, apg, amg
    real(dp), dimension(nlev) ::  w, g, f
    real(dp), dimension(nlev) :: tau_Ve_s
    real(dp), dimension(nlay) :: tau
    real(dp), dimension(nlev) :: tau_s, w_s, f_s, g_s
    real(dp), dimension(nlev) :: lam, u, N, gam, alp
    real(dp), dimension(nlev) :: R_b, T_b, R, T
    real(dp), dimension(nlev) :: Tf

    ! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    w(nlev) = w_surf
    g(nlev) = 0.0_dp

    ! Backscattering approximation
    f(:) = g(:)**2

    !! Do optical depth rescaling
    tau_Ve_s(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k) = tau_Ve(k+1) - tau_Ve(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_Ve_s(k+1) = tau_Ve_s(k) + tau_s(k)
    end do

    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_z**2)/(1.0_dp - lam(k)**2*mu_z**2)
      alp(k) = 0.75_dp * w_s(k) * mu_z * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_z**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_Ve_s(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_Ve_s(k)/mu_z,99.0_dp)
      Tf(k) = exp(-arg)

      apg = alp(k) + gam(k)
      amg = alp(k) - gam(k)

      R(k) = amg*(T_b(k)*Tf(k) - 1.0_dp) + apg*R_b(k)

      T(k) = apg*T_b(k) + (amg*R_b(k) - (apg - 1.0_dp))*Tf(k)

      R(k) = max(R(k), 0.0_dp)
      T(k) = max(T(k), 0.0_dp)
      R_b(k) = max(R_b(k), 0.0_dp)
      T_b(k) = max(T_b(k), 0.0_dp)

    end do

    !! Calculate downward flux
    do k = 1, nlay
      sw_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R(k+1)*R_b(k))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_down(nlev) = Tf(nlev)

    !! Calculate upward flux
    do k = 1, nlay
      sw_up(k) = (Tf(k)*R(k+1) + (T(k) - Tf(k))*R_b(k+1))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_up(nlev) = sw_down(nlev) * w_surf

    !! Scale with the incident flux
    sw_down(:) = sw_down(:) * mu_z * Finc
    sw_up(:) = sw_up(:) * mu_z * Finc

  end subroutine sw_grey_updown_adding

end module ts_isothermal_2_mod
