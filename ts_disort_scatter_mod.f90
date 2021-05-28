!!!
! Elspeth KH Lee - May 2021
! Two-stream DISORT version, modified by Xianyu Tan to include an internal heat source.
! Pros: Stable, performs accurate scattering calculations tried and tested, reliable model.
! Cons: Slower than other methods.
!!!

module ts_disort_scatter_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  public :: ts_disort_scatter
  private :: linear_log_interp

contains

  subroutine ts_disort_scatter(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp), dimension(nlev) :: Te

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc, utau
    real(dp), dimension(maxcly) :: gg, ssalb
    real(dp), dimension(maxulv) :: sw_net, lw_net
    real(dp) :: umu0, fbeam
    logical :: planck


    !! Find temperature at layer edges through linear interpolation and extrapolation
    do i = 2, nlay
      call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    end do
    Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    Te_0(0:nlay) = Te(1:nlev)

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      planck = .False.
      gg(1:nlay) = sw_g(:)
      ssalb(1:nlay) = sw_a(:)
      fbeam = F0
      umu0 = mu_z
      wvnmlo = 0.0_dp
      wvnmhi = 1.0e7_dp
      utau(1:nlev) = tau_Ve(:)
      do i = 1, nlay
        dtauc(i) = (tau_Ve(i+1) - tau_Ve(i))
      end do
      call CALL_TWOSTR (nlay,Te_0,gg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,sw_net)
    else
      sw_net(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    planck = .True.
    gg(1:nlay) = lw_g(:)
    ssalb(1:nlay) = lw_a(:)
    fbeam = 0.0_dp
    umu0 = 1.0_dp
    wvnmlo = 0.0_dp
    wvnmhi = 1.0e7_dp
    utau(1:nlev) = tau_IRe(:)
    do i = 1, nlay
      dtauc(i) = (tau_IRe(i+1) - tau_IRe(i))
    end do
    call CALL_TWOSTR (nlay,Te,gg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,lw_net)

    !! Net fluxes at each level
    net_F(:) = lw_net(1:nlev) + sw_net(1:nlev)

  end subroutine ts_disort_scatter

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

end module ts_disort_scatter_mod
