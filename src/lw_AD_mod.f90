!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Adding-Doubling method
!!!

module lw_AD_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: c1 = (hp * c_s) / kb
  real(dp), parameter :: c2 = c_s**2
  real(dp), parameter :: n2 = 2.0_dp * hp * c2

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  integer, parameter :: nmu = 1
  real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Gauss–Jacobi-5 quadrature for 2 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 2
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473988_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  !! Optical depth limit on optical depth of layers
  real(dp) :: dtau_lim = 0.01_dp

  private :: lw_doubling_adding, BB_integrate
  public :: lw_AD

contains


  subroutine lw_AD(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, a_surf, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: a_surf
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
    be(:) = sb * Te(:)**4
    be_int = sb * Tint**4

    !! Longwave flux calculation
    call lw_adding_doubling(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AD

  subroutine lw_adding_doubling(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
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
    real(dp), dimension(nlay) :: b1, b0, T
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


    ! TOA Reflection boundary
    R() = 
    ! TOA Transmission boundary
    T() = 
    ! TOA emissivity and slope
    y() =
    z() = 

    !! Begin loop from top of atmosphere to bottom to find transmision and emission and each layer boundary
    !! Perform doubling scheme as required when dtau > dtau_lim
    do k = 1, nlev

      if (dtau(k) < dtau_lim) then
        !! Safe to perform the adding method directly


      else
        !! Must perform the doubling method until dtau_lim is reached
        !! Its probably best to split into equal dtau sub-layers


      end if
    end do


  end subroutine lw_adding_doubling

  subroutine BB_integrate(n_b, Te, wn_e, be)
    implicit none

    integer, intent(in) :: n_b
    real(dp), intent(in) :: Te
    real(dp), dimension(n_b+1), intent(in) :: wn_e

    real(dp), dimension(n_b), intent(out) :: be

    integer :: ww, j, intitera
    real(dp), dimension(n_b+1) :: iB
    real(dp) :: x, x2, x3, itera, summ, dn

    !! Code for integrating the blckbody function between two wavenumbers
    !! This is a method that uses a sum convergence
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

      if (Te < 1e-6_dp) then
        be(:) = 0.0_dp
        return
      end if

      do ww = 1, n_b+1

        x = c1 * 100.0_dp * wn_e(ww)/ Te
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera > 150) then
          itera = 150
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,kind=dp)
          summ = summ +  exp(-min(real(j,kind=dp)*x,300.0_dp))* &
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        iB(ww) = n2 * (Te/c1)**(4) * summ
      end do

      do ww = 1, n_b
        be(ww) = max(iB(ww+1) - iB(ww),0.0_dp)
      end do

  end subroutine BB_integrate

end module lw_AD_mod