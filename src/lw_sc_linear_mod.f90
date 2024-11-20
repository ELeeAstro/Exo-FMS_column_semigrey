!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Oct 2021 : adding method & Bezier interpolation
!                - Aug 2023 : Change quadrature following Li (2000)
!                - Jun 2024 : Change quadrature following Hogan (2004)
! lw: Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
!     Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
!     Pros: Very fast, accurate at high optical depths, very stable
!     Cons: No lw scattering (but could combine with AA if needed)
!!!

module lw_sc_linear_mod
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

  private :: lw_shortchar_linear
  public :: lw_sc_linear

contains


  subroutine lw_sc_linear(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_shortchar_linear(nlay, nlev, be(:), be_int, tau_e(:),  ssa(:), gg(:), & 
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_sc_linear

  subroutine lw_shortchar_linear(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel, w0, hg, eps, dtau_a
    real(dp), dimension(nlay) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm, Gp, Bp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nlay) :: fc, sigma_sq, pmom2, c
    integer, parameter :: nstr = nmu*2

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)


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

      del(:) = dtau_a(:)/uarr(m)
      edel(:) = exp(-del(:))
      e0i(:) = 1.0_dp - edel(:)

      !! Prepare loop
      ! Olson & Kunasz (1987) linear interpolant parameters
      where (edel(:) > 0.999_dp)
        ! If we are in very low optical depth regime, then use an isothermal approximation
        Am(:) = (0.5_dp*(be(2:) + be(1:nlay)) * e0i(:))/be(1:nlay)
        Bm(:) = 0.0_dp
        Gp(:) = 0.0_dp
        Bp(:) = Am(:)
      elsewhere
        ! Use linear interpolants
        e1i(:) = del(:) - e0i(:)
        e1i_del(:) = e1i(:)/del(:) ! The equivalent to the linear in tau term

        Am(:) = e0i(:) - e1i_del(:) ! Am(k) = Gp(k), just indexed differently
        Bm(:) = e1i_del(:) ! Bm(k) = Bp(k), just indexed differently
        Gp(:) = Am(:)
        Bp(:) = Bm(:)
      end where

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k) + Bm(k)*be(k+1) ! TS intensity
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here we use the same condition but use intensity units to be consistent
      lw_up_g(nlev) = lw_down_g(nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bp(k)*be(k) + Gp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_shortchar_linear

end module lw_sc_linear_mod