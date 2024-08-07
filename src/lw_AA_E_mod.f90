!!!
! Elspeth KH Lee - Aug 2023 : Initial version
!
! lw: (Extended) Absorption Approximation following Li (2002) (alpha-nEAA) - This is the exoponential in tau method
!     Pros: Very fast method with LW scattering approximation, no matrix inversions
!     Cons: Not technically multiple scattering (can be quite innacurate at even moderate albedo)
!!!

module lw_AA_E_mod
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

  private :: lw_AA_exp
  public :: lw_AA_E

contains

  subroutine lw_AA_E(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_AA_exp(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AA_E

  subroutine lw_AA_exp(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, w0, hg, eps, dtau_a
    real(dp), dimension(nlay) :: bln, T, nup, nun
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nlay) :: fc, sigma_sq, pmom2, c
    integer, parameter :: nstr = nmu*2

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling
    where (g_in(:) >= 1e-6_dp)
      fc(:) = g_in(:)**(nstr)
      pmom2(:) = g_in(:)**(nstr+1)
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
      dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = w_in(:)
    end where

    hg(:) = g_in(:)

    !! Log B with tau function
    !where (dtau(:) < 1.0e-9_dp)
    !  bln(:) = 0.0_dp
    !elsewhere
      bln(:) = log(be(2:nlev)/be(1:nlay))/dtau(:)
    !end where

    !! modified co-albedo epsilon
    eps(:) = sqrt((1.0_dp - w0(:))*(1.0_dp - hg(:)*w0(:)))
    !eps(:) = (1.0_dp - w0(:))

    !! Absorption/modified optical depth for transmission function
    dtau_a(:) = eps(:)*dtau(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Efficency variables
      T(:) = exp(-dtau_a(:)/uarr(m))
      nup(:) = 1.0_dp/(1.0_dp + uarr(m)*bln(:)/eps(:))
      nun(:) = 1.0_dp/(1.0_dp - uarr(m)*bln(:)/eps(:))

      !! Begin two-stream loops
      !! Perform downward loop first - also calculate efficency variables
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*T(k) + &
          & nup(k) * (be(k+1) - be(k)*T(k))
      end do

      !! Perform upward loop
      !if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
      !  lw_up_g(nlev) = lw_down_g(nlev)*lw_a_surf + be_int
      !else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(nlev) = lw_down_g(nlev) + be_int
      !end if

      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*T(k) + &
          & nun(k) * (be(k) - be(k+1)*T(k))
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)
  

  end subroutine lw_AA_exp

end module lw_AA_E_mod