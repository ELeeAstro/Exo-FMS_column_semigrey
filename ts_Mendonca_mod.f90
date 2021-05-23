!!!
! Elspeth KH Lee - May 2021
! Two-stream method following Mendonca et al. methods used for Venus
! IN DEVELOPMENT - USE WITH UTMOST CAUTION, YOU WILL GET WRONG ANSWERS
!
!!!

module ts_Mendonca_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  public :: ts_Mendonca
  private :: lw_grey_updown, sw_grey_down, linear_log_interp

contains

  subroutine ts_Mendonca(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
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
    be(:) = (sb * Te(:)**4)  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4) ! Integrated planck function intensity for internal temperature
    call lw_grey_updown(nlay, nlev, be, be_int, tau_IRe(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

  end subroutine ts_Mendonca

  subroutine lw_grey_updown(nlay, nlev, be, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k
    real(dp), dimension(nlay) :: dtau, E, Es, r, mu, Tf

    !! Calculate dtau in each layer and prepare loop
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
      r(k) = 1.5_dp + 0.5_dp/(1.0_dp + 4.0_dp*dtau(k) + 10.0_dp*dtau(k)**2)
      mu(k) = 1.0_dp/r(k)
      Tf(k) = exp(-dtau(k)/mu(k))
      E(k) = ((be(k+1) - be(k) * Tf(k)) * dtau(k)) / (dtau(k) - mu(k)*log(be(k)/be(k+1)))
      Es(k) = ((be(k) - be(k+1) * Tf(k)) * dtau(k)) / (dtau(k) - mu(k)*log(be(k+1)/be(k)))
    end do

    !! Begin two-stream loops
    !! Perform downward loop first
    ! Top boundary condition - 0 flux downward from top boundary
    lw_down(1) = 0.0_dp
    do k = 1, nlay
      lw_down(k+1) = lw_down(k)*Tf(k) + Es(k)
    end do

    !! Perform upward loop
    ! Lower boundary condition - internal heat definition Fint = F_up - F_down
    lw_up(nlev) = lw_down(nlev) + be_int
    do k = nlay, 1, -1
      lw_up(k) = lw_up(k+1)*Tf(k) + E(k)
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

end module ts_Mendonca_mod