!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Feutrier
!     Pros: 
!     Cons: 
!!!

module lw_AA_L_mod
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
  

  private :: lw_Feutrier_method
  public :: lw_Feutrier

contains


  subroutine lw_Feutrier(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_Feutrier_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_Feutrier

  subroutine lw_Feutrier_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables - we follow the Hubeny indexing conventions
    integer :: d, l, k 
    real(dp), dimension(nlay) :: dtau, w0, hg
    integer, parameter :: n = 3

    real(dp), dimension(nlev,n,n) :: A, B, C
    real(dp), dimension(nlev,n) :: R, j
    real(dp), dimension(nlev) :: D, E

    w0(:) = w_in(:)
    hg(:) = g_in(:)

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    A(:,:,:) = 0.0_dp
    B(:,:,:) = 0.0_dp
    C(:,:,:) = 0.0_dp
    R(:) = 0.0_dp

    do d = 1, nlev


      if (d == 1) then

        ! Special boundary conditions
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k)
          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if


            A(d,l,k) = 0.0_dp
            C(d,l,k) = 2.0_dp * (mu(l)/dtau_32)**2 * kdel
            B(d,l,k) = (1.0_dp + (2.0_dp*mu(l)/dtau_32) + 2.0_dp*(mu(l)/dtau_32)**2) * kdel - alph(d,l) * phi(d,k)

          end do
        end do

      else if (d == nlev) then

        ! Special boundary conditions
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k) + (2.0_dp*mu(l)/dtau_mh) * be_int
          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if


            A(d,l,k) = 2.0_dp * (mu(l)/dtau_mh)**2 * kdel
            C(d,l,k) = 0.0_dp
            B(d,l,k) = (1.0_dp + (2.0_dp*mu(l)/dtau_mh) + 2.0_dp*(mu(l)/dtau_mh)**2) * kdel - alph(d,l) * phi(d,k)

          end do
        end do

      else
        !! Find the values of the large tridiagonal matrix
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k)

          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if

            A(d,l,k) = mu(l) / (dtau(d-1)*dtau(d)) * kdel
            C(d,l,k) = mu(l) / (dtau(d+1)*dtau(d)) * kdel
            B(d,l,k) = kdel +  A(d,l,k) + C(d,l,k) - alph(d,l) * w(d,k)

          end do
        end do 
      end if

    end do


    do l = 1, n
      ! Find D and E for this mu and solve for j
      do d = 1, nlev

      end do
    end do



  end subroutine lw_Feutrier_method

end module lw_AA_L_mod