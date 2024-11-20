!!!
! Elspeth KH Lee - July 2024 : Initial version
! lw: 
!!!

module lw_DFE_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node - 2 stream (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes - 4 stream (Hogan 2024)
  !integer, parameter :: nmu = 2
  !real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473981_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  !! Gauss–Jacobi-5 quadrature for 3 nodes - 6 stream (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1024922169_dp, 0.4417960320_dp, 0.8633751621_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0437820218_dp, 0.3875796738_dp, 0.5686383044_dp/)

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)

  private :: lw_discontinuous_finite_element
  public :: lw_DFE

contains


  subroutine lw_DFE(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    integer :: k
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nlev) :: be
    real(dp) :: be_int

    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_discontinuous_finite_element(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

    ! do k = 1, nlev
    !   print*, k, lw_net(k), lw_up(k), lw_down(k)
    ! end do

    ! stop

  end subroutine lw_DFE

  subroutine lw_discontinuous_finite_element(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), intent(in) :: be_int
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlay) :: del, a, b
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g, Ip, Im

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)/uarr(m)

      !! Find a and b for downward loop
      b(:) = del(:)*(del(:) + 1.0_dp)
      a(:) = del(:)**2 + 2.0_dp*del(:) + 2.0_dp

      !! Perform the Ip and Im calculation for the downward stream
      Im(1) = 0.0_dp ! No external radiation
      do k = 1, nlay
        Im(k+1) = (2.0_dp * Im(k) + del(k)*be(k) + b(k)*be(k+1))/a(k)
        Ip(k) = (2.0_dp*(del(k) + 1.0_dp)*Im(k) + b(k)*be(k) - del(k)*be(k+1))/a(k)
      end do
      Ip(nlev) = 0.0_dp

      lw_down_g(1) = Im(1)
      do k = 2, nlay
        lw_down_g(k) = (Im(k)*del(k) + Ip(k)*del(k-1))/(del(k) + del(k-1))
      end do
      lw_down_g(nlev) = Im(nlev)

      !! Perform the Ip and Im calculation for the upward stream
      Im(nlev) = be(nlev) + uarr(m)*(be(nlev) - be(nlev-1))/del(nlay)! Internal radiation
      do k = nlev, 2, -1
        Im(k-1) = (2.0_dp * Im(k) + del(k-1)*be(k) + b(k-1)*be(k-1))/a(k-1)
        Ip(k) = (2.0_dp*(del(k-1) + 1.0_dp)*Im(k) + b(k-1)*be(k) - del(k-1)*be(k-1))/a(k-1)
      end do
      Ip(1) = 0.0_dp

      lw_up_g(1) = Im(1)
      do k = 2, nlay
        lw_up_g(k) = (Im(k)*del(k) + Ip(k)*del(k-1))/(del(k) + del(k-1))
      end do
      lw_up_g(nlev) = Im(nlev)
      
      
      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_discontinuous_finite_element

end module lw_DFE_mod