!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Oct 2021 : adding method & Bezier interpolation
!                - Aug 2023 Change quadrature following Li (2000)
! sw: Adding layer method with scattering
! lw: Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
!     Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
!     Pros: Very fast, accurate at high optical depths, very stable
!     Cons: No lw scattering
!!!

module ts_short_char_mod_linear
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

  public :: ts_short_char_linear
  private :: lw_grey_updown_linear, sw_grey_updown_adding, linear_interp, Bezier_interp

contains

  subroutine ts_short_char_linear(Bezier, surf, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, &
    & mu_z, F0, Tint, AB, sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr, net_Fs)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier, surf
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Ts, F0, Tint, AB, sw_a_surf, lw_a_surf
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

    !! Output variables
    real(dp), intent(out) :: olr, asr, net_Fs
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: lpe
    real(dp), dimension(nlay) :: lTl, lpl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    ! Log the layer values and pressure edges for more accurate interpolation
    lTl(:) = log10(Tl(:))
    lpl(:) = log10(pl(:))
    lpe(:) = log10(pe(:))

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(lpl(i-1:i+1), lTl(i-1:i+1), 3, lpe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(lpl(nlay-2:nlay), lTl(nlay-2:nlay), 3, lpe(nlay), Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_interp(lpe(i), lpl(i-1), lpl(i), lTl(i-1), lTl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = (log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/10.0_dp**Te(2)))
    Te(nlev) = (log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/10.0_dp**Te(nlay)))
    ! De-log the temperature levels (edges)
    Te(:) = 10.0_dp**Te(:)

    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve(:), mu_z(:), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    if (surf .eqv. .True.) then
      be_int = (sb * Ts**4)/pi
    else
      be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    end if

    call lw_grey_updown_linear(surf, nlay, nlev, be, be_int, lw_a, lw_g, lw_a_surf, tau_IRe(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Net surface flux (for surface temperature evolution)
    !! We have to define positive as downward (heating) and cooling (upward) in this case
    net_Fs = sw_down(nlev) + lw_down(nlev) - lw_up(nlev)

    !! Output the olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_short_char_linear

  subroutine lw_grey_updown_linear(surf, nlay, nlev, be, be_int, lw_a, lw_g, lw_a_surf, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), dimension(nlay), intent(in) :: lw_a, lw_g
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp), dimension(nlay) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm, Gp, Bp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:) - tau_IRe(1:nlay)

    !! Absorption approxmation
    dtau(:) =  dtau(:) * (1.0_dp - lw_a(:))

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)/uarr(m)
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
      if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        lw_up_g(nlev) = lw_down_g(nlev)*lw_a_surf + be_int
      else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(nlev) = lw_down_g(nlev) + be_int
      end if

      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bp(k)*be(k) + Gp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * w(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi (in this GJ weighting scheme)
    lw_down(:) = pi * lw_down(:)
    lw_up(:) = pi * lw_up(:)

  end subroutine lw_grey_updown_linear

  subroutine sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve, mu_z
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables
    integer :: k
    real(dp) :: lamtau, e_lamtau, arg, apg, amg
    real(dp), dimension(nlev) ::  w, g, f
    real(dp), dimension(nlev) :: tau_Ve_s
    real(dp), dimension(nlay) :: tau
    real(dp), dimension(nlev) :: tau_s, w_s, g_s
    real(dp), dimension(nlev) :: lam, u, N, gam, alp
    real(dp), dimension(nlev) :: R_b, T_b, R, T
    real(dp), dimension(nlev) :: Tf
    real(dp), dimension(nlev) :: cum_trans

    ! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then

      if (mu_z(nlev) == mu_z(1)) then
        ! No zenith correction, use regular method
        sw_down(:) = Finc * mu_z(nlev) * exp(-tau_Ve(:)/mu_z(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_Ve(1)/mu_z(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_Ve(k+1) - tau_Ve(k))/mu_z(k+1)
        end do
        do k = 1, nlev
          sw_down(k) = Finc * mu_z(nlev) * exp(-cum_trans(k))
        end do
      end if

      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo

      return

    end if

    w(nlev) = w_surf
    g(nlev) = 0.0_dp

    ! Backscattering approximation
    f(:) = g(:)**2

    !! Do optical depth rescaling
    tau_Ve_s(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k) = tau_Ve(k+1) - tau_Ve(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_Ve_s(k+1) = tau_Ve_s(k) + tau_s(k)
    end do

    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_z(k)**2)/(1.0_dp - lam(k)**2*mu_z(k)**2)
      alp(k) = 0.75_dp * w_s(k) * mu_z(k) * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_z(k)**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_Ve_s(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_Ve_s(k)/mu_z(k),99.0_dp)
      Tf(k) = exp(-arg)

      apg = alp(k) + gam(k)
      amg = alp(k) - gam(k)

      R(k) = amg*(T_b(k)*Tf(k) - 1.0_dp) + apg*R_b(k)

      T(k) = apg*T_b(k) + (amg*R_b(k) - (apg - 1.0_dp))*Tf(k)

      R(k) = max(R(k), 0.0_dp)
      T(k) = max(T(k), 0.0_dp)
      R_b(k) = max(R_b(k), 0.0_dp)
      T_b(k) = max(T_b(k), 0.0_dp)

    end do

    !! Calculate downward flux
    do k = 1, nlay
      sw_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R(k+1)*R_b(k))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_down(nlev) = Tf(nlev)

    !! Calculate upward flux
    do k = 1, nlay
      sw_up(k) = (Tf(k)*R(k+1) + (T(k) - Tf(k))*R_b(k+1))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_up(nlev) = sw_down(nlev) * w_surf

    !! Scale with the incident flux
    sw_down(:) = sw_down(:) * mu_z(nlev) * Finc
    sw_up(:) = sw_up(:) * mu_z(nlev) * Finc

  end subroutine sw_grey_updown_adding

  ! Perform linear interpolation
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  subroutine Bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine Bezier_interp

end module ts_short_char_mod_linear
