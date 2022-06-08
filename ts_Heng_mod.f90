!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Dec 2021 : adding method & Bezier interpolation
! sw: Adding layer method with scattering
! lw: Two-stream method following Heng et al. papers
!     We follow the Malik et al. (2017) method using sub-layers to calculate midpoint fluxes
!     Pros: Easy to understand and convert from theory, no mu integration, variable diffusive factor
!     Cons: Slower than other methods (uses sub-layers and eone calculations), no scattering
!     NOTE: Various ways to calculate the diffusion factor as function of optical depth
!!!

module ts_Heng_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  real(dp), parameter :: D = 1.66_dp  ! Diffusion factor
  real(dp), parameter :: eps = 1.0_dp/D  ! Diffusion factor

  real(dp), dimension(0:5) :: Aj = (/-0.57721566_dp, 0.99999193_dp, -0.24991055_dp, 0.05519968_dp, 0.00976004_dp, 0.00107857_dp /)
  real(dp), dimension(0:4) :: Bj = (/1.0_dp, 8.5733287401_dp, 18.059016973_dp, 8.6347608925_dp, 0.2677737343_dp /)
  real(dp), dimension(0:4) :: Cj = (/1.0_dp, 9.5733223454_dp, 25.6329561486_dp, 21.0996530827_dp, 3.9584969228_dp /)

  public :: ts_Heng
  private :: lw_grey_updown, sw_grey_updown_adding, linear_log_interp, bezier_interp, e1_ap

contains

  subroutine ts_Heng(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, tau_IRl, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, sw_a_surf, net_F, olr, asr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB, sw_a_surf
    real(dp), dimension(nlay), intent(in) :: Tl, pl, tau_IRl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g

    !! Output variables
    real(dp),  intent(out) :: olr, asr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: bl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(pl(i-1:i+1), Tl(i-1:i+1), 3, pe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(pl(nlay-2:nlay), Tl(nlay-2:nlay), 3, pe(nlay), Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve(:), mu_z, sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = sb * Te(:)**4/pi  ! Integrated planck function intensity at levels
    bl(:) = sb * Tl(:)**4/pi  ! Integrated planck function intensity at layers
    be_int = sb * Tint**4/pi ! Integrated planck function intensity for internal temperature
    call lw_grey_updown(nlay, nlev, be, bl,  be_int, tau_IRe(:), tau_IRl(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_Heng

  subroutine lw_grey_updown(nlay, nlev, be, bl, be_int, tau_IRe, tau_IRl, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), dimension(nlay), intent(in) :: bl, tau_IRl
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k
    real(dp) :: Fmid, Tp, Bp, dtau
    !real(dp), dimension(nlay) :: dtau, Tp, Bp

    real(dp), external :: eone

    lw_down(1) = 0.0_dp

   do k = 1, nlay

      ! From upper to mid
      !! delta tau
      dtau = (tau_IRl(k) - tau_IRe(k))
      !! Transmission function
      !Tp = exp(-D*dtau)
      Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*eone(dtau)
      !Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*e1_ap(dtau)
      !! Linear in tau approximation
      if (dtau < 1e-6_dp) then
        Bp = 0.0_dp
      else
        Bp = (bl(k) - be(k))/dtau
      end if

      !! Flux expression
      !Fmid = lw_down_b(k)*Tp + pi*be(k)*(1.0_dp - Tp) ! Isothermal approximation
      Fmid =  lw_down(k)*Tp + pi*be(k)*(1.0_dp - Tp) + &
       & pi*Bp * (-2.0_dp/3.0_dp*(1.0_dp - exp(-dtau)) + dtau*(1.0_dp - Tp/3.0_dp))

      ! From mid to lower
      dtau = (tau_IRe(k+1) - tau_IRl(k))
      !Tp = exp(-D*dtau)
      Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*eone(dtau)
      !Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*e1_ap(dtau)
      if (dtau < 1e-6_dp) then
        Bp = 0.0_dp
      else
        Bp = (be(k+1) - bl(k))/dtau
      end if

      !lw_down_b(k+1) =  Fmid*Tp + pi*bl(k)*(1.0_dp - Tp) ! Isothermal approximation
      lw_down(k+1) =  Fmid*Tp + pi*bl(k)*(1.0_dp - Tp) + &
      & pi*Bp * (-2.0_dp/3.0_dp*(1.0_dp - exp(-dtau)) + dtau*(1.0_dp - Tp/3.0_dp))

    end do


    ! Upward boundary condition - NOTE: contains intenal flux contribution
    lw_up(nlev) = lw_down(nlev) + pi*be_int !pi*(bsurf)

    !! Peform upward loop
    do k = nlay, 1, -1

      ! From lower to mid
      dtau = (tau_IRe(k+1) - tau_IRl(k))
      !Tp = exp(-D*dtau)
      Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*eone(dtau)
      !Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*e1_ap(dtau)
      if (dtau < 1e-6_dp) then
        Bp = 0.0_dp
      else
        Bp = (be(k+1) - bl(k))/dtau
      end if

      !Fmid = lw_up_b(k+1)*Tp + pi*be(k+1)*(1.0_dp - Tp) ! Isothermal approximation
      Fmid = lw_up(k+1)*Tp + pi*be(k+1)*(1.0_dp - Tp) + &
      & pi*Bp * (2.0_dp/3.0_dp*(1.0_dp - exp(-dtau)) - dtau*(1.0_dp - Tp/3.0_dp))

      ! From mid to upper
      dtau = (tau_IRl(k) - tau_IRe(k))
      !Tp = exp(-D*dtau)
      Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*eone(dtau)
      !Tp = (1.0_dp - dtau)*exp(-dtau) + dtau**2*e1_ap(dtau)
      if (dtau < 1e-6_dp) then
        Bp = 0.0_dp
      else
        Bp = (bl(k) - be(k))/dtau
      end if

      !lw_up_b(k) = Fmid*Tp + pi*bl(k)*(1.0_dp - Tp) ! Isothermal approximation
      lw_up(k) = Fmid*Tp + pi*bl(k)*(1.0_dp - Tp) + &
      & pi*Bp * (2.0_dp/3.0_dp*(1.0_dp - exp(-dtau)) - dtau*(1.0_dp - Tp/3.0_dp))

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

  function e1_ap(dtau) result(e1)
    implicit none

    real(dp), intent(in) :: dtau
    real(dp) :: e1

    integer :: j
    real(dp) :: Asum, Bsum, Csum
    real(dp) :: dtau4

    if (dtau < 1.0_dp) then
      Asum = 0.0_dp
      do j = 0, 5
        Asum = Asum + Aj(j)*dtau**j
      end do
      e1 = -log(dtau) + Asum
    else
      Bsum = 0.0_dp
      Csum = 0.0_dp
      do j = 0, 4
        dtau4 = dtau**(4-j)
        Bsum = Bsum + Bj(j)*dtau4
        Csum = Csum + Cj(j)*dtau4
      end do
      e1 = (exp(-dtau)/dtau) * (Bsum/Csum)
    end if

  end function e1_ap

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w < min(wlim,wlim1) .or. w > max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w < min(wlim,wlim1) .or. w > max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_Heng_mod
