!!!
! Elspeth KH Lee - March 2024 - special max efficency version
!!!

module ts_maxeff_mod
  use, intrinsic :: iso_fortran_env
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
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes (Hogan 2024)
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473988_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  public :: ts_maxeff
  private :: lw_maxeff, sw_maxeff, linear_log_interp, bezier_interp

contains

  subroutine ts_maxeff(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & net_F, olr, asr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, Tint, AB, mu_z
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe

    !! Output variables
    real(dp), intent(out) :: olr, asr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: lpe
    real(dp), dimension(nlay) :: lTl, lpl
    real(dp), dimension(nlev) :: sw_down, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then

      ! Log the layer values and pressure edges for more accurate interpolation
      lTl(:) = log10(Tl(:))
      lpl(:) = log10(pl(:))
      lpe(:) = log10(pe(:))

      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(lpl(i-1:i+1), lTl(i-1:i+1), 3, lpe(i), Te(i))
        Te(i) = 10.0_dp**(Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(lpl(nlay-2:nlay), lTl(nlay-2:nlay), 3, lpe(nlay), Te(nlay))
      Te(nlay) = 10.0_dp**(Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_maxeff(nlev, Finc, tau_Ve(:), mu_z, sw_down(:))
    else
      sw_down(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    call lw_maxeff(nlay, nlev, be, be_int, tau_IRe(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = 0.0_dp - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1)

  end subroutine ts_maxeff

  subroutine lw_maxeff(nlay, nlev, be, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlay) :: b1, b0, T
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    dtau(:) = tau_IRe(2:nlev) - tau_IRe(1:nlay)

    !! Linear B with tau function
    where (dtau(:) <= 1.0e-6_dp)
      b1(:) = 0.0_dp
      b0(:) = 0.5_dp*(be(2:nlev) + be(1:nlay))
    elsewhere
      b1(:) = (be(2:nlev) - be(1:nlay))/dtau(:) 
      b0(:) = be(1:nlay)
    end where

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Transmission function
      T(:) = exp(-dtau(:)/uarr(m))

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*T(k) + b0(k)*(1.0_dp - T(k)) + b1(k)*(uarr(m)*T(k)+dtau(k)-uarr(m)) ! TS intensity
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_up - F_down
      lw_up_g(nlev) = lw_down_g(nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*T(k) + b0(k)*(1.0_dp - T(k)) + b1(k)*(uarr(m)-(dtau(k)+uarr(m))*T(k)) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * w(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    lw_down(:) = pi * lw_down(:)
    lw_up(:) = pi * lw_up(:)

  end subroutine lw_maxeff

  subroutine sw_maxeff(nlev, Finc, tau_Ve, mu_z, sw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_Ve

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down

    !! Direct beam only
    sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)

  end subroutine sw_maxeff


  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: ly1, ly2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, wh, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if ((x > xi(1)) .and. (x < xi(2))) then
      ! left hand side interpolation
      !print*,'left'
      wh = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
        wh = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (wh*dy/dx + (1.0_dp - wh)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      wh = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
        wh = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (wh*dy1/dx1 + (1.0_dp - wh)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_maxeff_mod
