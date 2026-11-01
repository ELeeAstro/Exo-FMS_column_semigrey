!!!
! Elspeth KH Lee - Jun 2024 : Overhauled
!
! sw: Adding method for shortwave radiation, used in various Earth GCMs (e.g. Mendonça et al. 2015)
!     Pros: Very fast method with approximate scattering, decently accurate
!     Cons: Not multiple scattering
!!!


module sw_adding_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_approx_adding
  public :: sw_adding

contains

  subroutine sw_adding(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev) :: mu_z
    real(dp) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Shortwave flux calculation
    if ((mu_z(nlev) > 0.0_dp) .and. (Finc > 0.0_dp)) then
      call sw_approx_adding(nlay, nlev, Finc, mu_z(:), tau_e(:), ssa(:), gg(:), &
        & sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_adding


  subroutine sw_approx_adding(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    integer :: k
    real(dp) :: lamtau, e_lamtau, arg, apg, amg
    real(dp) :: tau_lo, N_lo, Tf_lo
    real(dp), dimension(nlev) ::  w, g, f
    real(dp), dimension(nlev) :: tau_e
    real(dp), dimension(nlay) :: tau
    real(dp), dimension(nlev) :: tau_s, w_s, g_s
    real(dp), dimension(nlev) :: lam, u, N, gam, alp
    real(dp), dimension(nlev) :: R_b, T_b, T
    real(dp), dimension(nlev) :: R_lo, Rb_lo, Tb_lo
    real(dp), dimension(nlev) :: Tf
    real(dp), dimension(nlev) :: cum_trans

    !! Add a black gas-giant lower-boundary entry.
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        do k = 1, nlev
          flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
        end do
      end if

      flx_up(:) = 0.0_dp

      return

    end if

    !! Backscattering approximation
    f(:) = g(:)**2

    !! Do optical depth rescaling
    tau_e(1) = tau_in(1)
    do k = 1, nlay
      tau(k) = tau_in(k+1) - tau_in(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_e(k+1) = tau_e(k) + tau_s(k)
    end do

    !! Perform approximate adding method.
    !! At every level we need the properties of BOTH slabs meeting there: the
    !! slab above (optical depth tau_e(k), measured down from the top) and the
    !! slab below (the remaining optical depth tau_e(nlev) - tau_e(k)). The two
    !! are then combined with the adding formulae further down.
    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_in(k)**2)/(1.0_dp - lam(k)**2*mu_in(k)**2)
      alp(k) = 0.75_dp * w_s(k) * mu_in(k) * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_in(k)**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      apg = alp(k) + gam(k)
      amg = alp(k) - gam(k)

      !! Slab above level k - diffuse reflectance R_b, diffuse transmittance
      !! T_b, direct transmittance Tf and total (direct + diffuse) transmittance T
      lamtau = min(lam(k)*tau_e(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_e(k)/mu_in(k),99.0_dp)
      Tf(k) = exp(-arg)

      T(k) = apg*T_b(k) + (amg*R_b(k) - (apg - 1.0_dp))*Tf(k)

      !! Slab below level k - direct-beam reflectance R_lo and diffuse
      !! reflectance Rb_lo. These are what the level sees looking downward,
      !! so they must be built from the optical depth that remains below it.
      tau_lo = tau_e(nlev) - tau_e(k)

      lamtau = min(lam(k)*tau_lo,99.0_dp)
      e_lamtau = exp(-lamtau)

      N_lo = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      Rb_lo(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N_lo
      Tb_lo(k) = 4.0_dp * u(k)/N_lo

      arg = min(tau_lo/mu_in(k),99.0_dp)
      Tf_lo = exp(-arg)

      R_lo(k) = amg*(Tb_lo(k)*Tf_lo - 1.0_dp) + apg*Rb_lo(k)

      T(k) = max(T(k), 0.0_dp)
      R_b(k) = max(R_b(k), 0.0_dp)
      T_b(k) = max(T_b(k), 0.0_dp)
      R_lo(k) = max(R_lo(k), 0.0_dp)
      Rb_lo(k) = max(Rb_lo(k), 0.0_dp)

    end do

    !! Add the slab above and the slab below at each level. At k = 1 the upper
    !! slab is empty (Tf = T = 1, R_b = 0) so flx_up(1) reduces to R_lo(1), the
    !! reflectance of the whole column - i.e. the planetary albedo. At k = nlev
    !! the lower slab is empty (R_lo = Rb_lo = 0) so flx_up(nlev) = 0, the
    !! gas-giant lower boundary condition.
    do k = 1, nlev
      flx_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R_lo(k)*R_b(k))/(1.0_dp - R_b(k)*Rb_lo(k))
      flx_up(k) = (Tf(k)*R_lo(k) + (T(k) - Tf(k))*Rb_lo(k))/(1.0_dp - R_b(k)*Rb_lo(k))
    end do

    !! Scale with the incident flux
    flx_down(:) = flx_down(:) * mu_in(nlev) * F0_in
    flx_up(:) = flx_up(:) * mu_in(nlev) * F0_in

  end subroutine sw_approx_adding

end module sw_adding_mod
