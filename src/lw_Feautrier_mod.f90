!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Feautrier
!     Pros: 
!     Cons: 
!!!

module lw_Feautrier_mod
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
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  integer, parameter :: nmu = 1
  real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 2
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473988_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)


  private :: lw_Feautrier_method
  public :: lw_Feautrier

contains


  subroutine lw_Feautrier(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_Feautrier_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_Feautrier

  subroutine lw_Feautrier_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
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
    real(dp), dimension(nlay) :: dtau, w0, hg

    real(dp), dimension(nlev,nmu) :: A, B, C
    real(dp), dimension(nlev) :: R

    w0(:) = w_in(:)
    hg(:) = g_in(:)

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    ! Only need the diagonal elements of the matrix, so arrays are nlev by nmu
    A(:,:) = 0.0_dp
    B(:,:) = 0.0_dp
    C(:,:) = 0.0_dp

    R(:) = 0.0_dp

    do k = 1, nlev


      if (k == 1) then

        ! Thermal source function at top of atmosphere
        R(k) = (1.0_dp - w0(k))*be(k)

        ! Special upper boundary conditions
        do m = 1, nmu

        end do

      else if (k == nlev) then

        ! Thermal source at lower boundary including internal temperature
        R(k) = (1.0_dp - w0(k))*be(k) + be_int

        ! Special lower boundary conditions
        do m = 1, nmu
        end do

      else

        ! Thermal source at level
        R(k) = (1.0_dp - w0(k))*be(k)
        do m = 1, nmu
        end do 

      end if


      ! Solve the tridiagonal system for this level

    end do

  end subroutine lw_Feautrier_method

end module lw_Feautrier_mod