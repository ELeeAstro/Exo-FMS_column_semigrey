!!!
! Elspeth KH Lee - Jun 2024
! lw: two-stream disort method with modifications from Xianyu Tan.
!     Pros: Accurate and reasonably fast multiple scattering method
!     Cons: Two-stream only, can be inaccurate
!!!

module lw_disort_ts_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private 
  public :: lw_disort_ts

contains


  subroutine lw_disort_ts(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr


    !! Work variables
    integer :: i
    real(dp), dimension(nlev) :: Te

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc
    real(dp), dimension(maxcly) :: ggg, ssalb
    real(dp), dimension(maxulv) :: lw_net_d, lw_up_d, lw_down_d, utau
    real(dp) :: umu0, fbeam
    logical :: planck


    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    Te_0(:) = 0.0_dp
    Te_0(0:nlay) = Te(1:nlev)

    !! Longwave two-stream flux calculation
    planck = .True.
    fbeam = 0.0_dp
    umu0 = 0.0_dp
    lw_net_d(:) = 0.0_dp
    lw_up_d(:) = 0.0_dp
    lw_down_d(:) = 0.0_dp
    ggg(:) = 0.0_dp
    ssalb(:) = 0.0_dp
    utau(:) = 0.0_dp
    dtauc(:) = 0.0_dp

    wvnmlo = 0.0_dp
    wvnmhi = 1e7_dp

    ggg(1:nlay) = gg(:)
    ssalb(1:nlay) = ssa(:)
    utau(1:nlev) = tau_e(:)
    do i = 1, nlay
      dtauc(i) = (tau_e(i+1) - tau_e(i))
    end do

    call call_twostr(nlay,Te_0,ggg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0, &
      & lw_net_d(:),lw_up_d(:),lw_down_d(:))

    lw_net(:) = lw_net_d(1:nlev)
    lw_up(:) = lw_up_d(1:nlev)
    lw_down(:) = lw_down_d(1:nlev)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_disort_ts

end module lw_disort_ts_mod