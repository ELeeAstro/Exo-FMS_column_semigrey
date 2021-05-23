!!!
! Elspeth KH Lee - May 2021
! Two-stream method following Heng et al. papers
! We follow the Malik et al. (2017) method using sub-layers to calculate midpoint fluxes
! Pros: Easy to understand and convert from theory, no mu integration, variable diffusive factor
! Cons: Slower than other methods (uses sub-layers and eone calculations), no scattering
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

  public :: ts_Heng
  private :: lw_grey_updown, sw_grey_down, linear_log_interp

contains

  subroutine ts_Heng(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, tau_IRl, mu_z, F0, Tint, AB, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl, tau_IRl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: bl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through linear interpolation and extrapolation
    do i = 2, nlay
      call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    end do
    Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_down(nlev, Finc, tau_Ve(:), mu_z, sw_down(:))
    else
      sw_down(:) = 0.0_dp
    end if
     sw_up(:) = 0.0_dp ! sw_up is zero since we don't have shortwave scattering in this mode

    !! Longwave two-stream flux calculation
    be(:) = sb * Te(:)**4/pi  ! Integrated planck function intensity at levels
    bl(:) = sb * Tl(:)**4/pi  ! Integrated planck function intensity at layers
    be_int = sb * Tint**4/pi ! Integrated planck function intensity for internal temperature
    call lw_grey_updown(nlay, nlev, be, bl,  be_int, tau_IRe(:), tau_IRl(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

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

  subroutine sw_grey_down(nlev, Finc, tau_V, mu_z, sw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_V

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_down

    sw_down(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

  end subroutine sw_grey_down

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

end module ts_Heng_mod
