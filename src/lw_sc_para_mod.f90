!!!
! Elspeth KH Lee - July 2024 : Initial version
! lw: 
!!!

module lw_sc_para_mod
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
  !! For HJ's we actually care about stratospheric heating rates, so Gauss–Laguerre quadrature is generally best for 4+ stream, 
  !! even for cloudy regions (Hogan 2024).
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Laguerre quadrature for 2 nodes (Hogan 2024)
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.1813898346_dp, 0.7461018061_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.1464466094_dp, 0.8535533906_dp/)

  !! Gauss–Laguerre quadrature for 3 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.0430681066_dp, 0.3175435896_dp, 0.8122985952_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0103892565_dp, 0.2785177336_dp, 0.7110930099_dp/)

  !! Gauss–Laguerre quadrature for 4 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.0091177205_dp, 0.1034869099_dp, 0.4177464746_dp, 0.8510589811_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0005392947_dp, 0.0388879085_dp, 0.3574186924_dp, 0.6031541043_dp /)

  private :: lw_shortchar_para
  public :: lw_sc_para

contains


  subroutine lw_sc_para(nlay, nlev, Tl, pl, pe, tau_e, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_e
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
    call lw_shortchar_para(nlay, nlev, be(:), be_int, tau_e(:), lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

    do k = 1, nlev
      print*, k, lw_net(k), lw_up(k), lw_down(k)
    end do

    stop

  end subroutine lw_sc_para

  subroutine lw_shortchar_para(nlay, nlev, be, be_int, tau_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp), dimension(nlay) :: del
    real(dp), dimension(nlay) :: Am, Bm, Cm, Ap, Bp, Cp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)/uarr(m)
      edel(:) = exp(-del(:))

      do k = 1, nlay
      end do


      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k-1) + Bm(k)*be(k+1)  + Cm(k)*be(k+1)! TS intensity
      end do

      lw_up_g(nlev) = lw_down_g(nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Ap(k)*be(k-1) + Bp(k)*be(k+1)  + Cp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_shortchar_para

end module lw_sc_para_mod