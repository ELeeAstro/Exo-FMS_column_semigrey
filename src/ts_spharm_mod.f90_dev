!!!
! Elspeth K.H. Lee - 
! sw: 
! lw: 
! Pros:
! Cons:
!!!

module ts_spharm_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Use Two-Term HG function for sw or lw
  logical, parameter :: TTHG_sw = .False.
  logical, parameter :: TTHG_lw = .False.

  public :: ts_sph
  private :: lw_sph, sw_sph, linear_log_interp, bezier_interp

contains

  subroutine ts_sph(surf, Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr, net_Fs)
    implicit none

    !! Input variables
    logical, intent(in) :: surf, Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, Tint, AB, sw_a_surf, lw_a_surf, Ts
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
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(pl(i-1:i+1), Tl(i-1:i+1), 3, pe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(pl(nlay-2:nlay), Tl(nlay-2:nlay), 3, pe(nlay), Te(nlay))
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
    if (mu_z(nlev) > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_sph(nlay, nlev, Finc, tau_Ve(:), mu_z(nlev), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
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

    call lw_sph(surf, nlay, nlev, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

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

  end subroutine ts_sph

  subroutine lw_sph(surf, nlay, nlev, be, be_int, tau_IRe, ww, gg, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf 
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), dimension(nlay), intent(in) :: ww, gg
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, i, j
    real(dp), dimension(nlay) :: sigma_sq, pmom2, c
    integer, parameter :: nstr = 4
    real(dp), dimension(nlay) :: w0, hg, fc
    real(dp), dimension(nlay) :: dtau

    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:) - tau_IRe(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling
    fc(:) = gg(:)**(nstr)
    pmom2(:) = gg(:)**(nstr+1)

    where (fc(:) /=  pmom2(:))
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = ww(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*ww(:)))
      dtau(:) = (1.0_dp - ww(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = ww(:)
      fc(:) = 0.0_dp
    end where

    hg(:) = gg(:)



  end subroutine lw_sph

  subroutine sw_sph(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables and arrays
    integer :: k
    real(dp), dimension(nlay) :: w0, dtau, hg
    real(dp), dimension(nlev) :: tau, dir
    real(dp) :: f0
    real(dp) :: om0, om1, om2, om3
    real(dp) :: a0, a1, a2, a3, b0, b1, b2, b3
    real(dp) :: e1, e2
    real(dp) :: beta, gam, lam1, lam2, R1, R2, Q1, Q2, S1, S2
    real(dp) :: eta0, eta1, eta2, eta3, del0, del1, del2, del3, delta
    real(dp), dimension(nlay) :: z1p, z1m, z2p, z2m
    real(dp), dimension(nlay) :: p1p, p1m, p2p, p2m
    real(dp), dimension(nlay) :: q1p, q1m, q2p, q2m
    real(dp), dimension(4,4,nlay) :: f
    real(dp), dimension(4,nlay) :: z_d, z_u

    real(dp), dimension(11,4*nlay) :: Mb ! Banded matrix
    real(dp), dimension(4*nlay) :: BB, F_bot 
    real(dp), dimension(4*nlev,4*nlay) :: FF
    real(dp), dimension(4*nlev) :: GG 
    real(dp) :: G_bot, b_top, b_surf

    real(dp), dimension(nlay) :: fc, sigma_sq, pmom2, c
    integer, parameter :: nstr = 4
    real(dp), parameter :: eps_20 = 1.0e-20_dp

    real(dp), dimension(nlay) :: dtr
    real(dp) :: hg2, alp

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-6_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_Ve(2:) - tau_Ve(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling
    fc(:) = g_in(:)**(nstr)
    pmom2(:) = g_in(:)**(nstr+1)

    where (fc(:) /=  pmom2(:))
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
      dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = w_in(:)
      fc(:) = 0.0_dp
    end where

    hg(:) = g_in(:)

    !! Reform edge optical depths
    tau(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Layer transmission
    dtr(:) = exp(-dtau(:)/mu_z)

    !! Direct beam transmission to layer
    dir(:) = exp(-tau(:)/mu_z)

    !! Start Spherical Harmonic sw method
    !! We follow the Rooney diagonal matrix method, 
    !! but also take inspiration from the SDA method of Zhang and Li (2013)
    do k = 1, nlay

      !! Inverse zenith angle
      f0 = 1.0_dp/mu_z

      !! Omega Legendre polynomial coefficents - scale with delta-M+
      if (hg(k) /= 0.0_dp) then
        if (TTHG_sw .eqv. .False.) then
          ! Use HG phase function
          om0 = 1.0_dp
          om1 = 3.0_dp * (hg(k) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * (hg(k)**2 - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * (hg(k)**3 - fc(k))/(1.0_dp - fc(k))
        else
          ! Use TTHG phase function with default parameters
          hg2 = hg(k)/2.0_dp
          alp = 1.0_dp - hg2**2
          om0 = 1.0_dp
          om1 = 3.0_dp * ((alp*hg(k) + (1.0_dp - alp)*hg2) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * ((alp*hg(k)**2 + (1.0_dp - alp)*hg2**2) - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * ((alp*hg(k)**3 + (1.0_dp - alp)*hg2**3) - fc(k))/(1.0_dp - fc(k))
        end if
      else
        ! Use Rayleigh scattering phase function for isotropic scattering
        om0 = 1.0_dp
        om1 = 0.0_dp
        om2 = 0.5_dp
        om3 = 0.0_dp
      end if

      !! Find the a coefficents
      a0 = 1.0_dp - w0(k)*om0 + eps_20
      a1 = 3.0_dp - w0(k)*om1 + eps_20
      a2 = 5.0_dp - w0(k)*om2 + eps_20
      a3 = 7.0_dp - w0(k)*om3 + eps_20

      !! Find the b coefficents - normalise Finc to 1 here
      b0 = w0(k)*om0 / fourpi
      b1 = w0(k)*om1 * -(mu_z) / fourpi
      b2 = 0.5_dp * w0(k)*om2 * (3.0_dp * mu_z**2 - 1.0_dp) / fourpi
      b3 = 0.5_dp * w0(k)*om3 * (5.0_dp * -mu_z**3 - 3.0_dp*-(mu_z)) / fourpi

      !! Find beta and gamma
      beta = a0*a1 + (4.0_dp/9.0_dp)*a0*a3 + (1.0_dp/9.0_dp)*a2*a3
      gam = (1.0_dp/9.0_dp)*a0*a1*a2*a3

      !! Find lambda values
      lam1 = sqrt((beta + sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp)
      lam2 = sqrt((beta - sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp)

      !! Find exponential values
      e1 = exp(-lam1*dtau(k))
      e2 = exp(-lam2*dtau(k))   

      !! Find S, R and Q coefficents
      S1 = -3.0_dp/2.0_dp * (a0*a1/lam1 - lam1)/a3
      S2 = -3.0_dp/2.0_dp * (a0*a1/lam2 - lam2)/a3
      R1 = -a0/lam1
      R2 = -a0/lam2
      Q1 = 0.5_dp * (a0*a1/lam1**2 - 1.0_dp)
      Q2 = 0.5_dp * (a0*a1/lam2**2 - 1.0_dp)

      !! Find the delta values
      delta = 9.0_dp*(f0**4 - beta*f0**2 + gam)
      del0 = (a1*b0 - b1*f0)*(a2*a3 - 9.0_dp*f0**2) + 2.0_dp*f0**2*(a3*b2 - 2.0_dp*a3*b0 - 3.0_dp*b3*f0)
      del1 = (a0*b1 - b0*f0)*(a2*a3 - 9.0_dp*f0**2) - 2.0_dp*a0*f0*(a3*b2 - 3.0_dp*b3*f0)
      del2 = (a3*b2 - 3.0_dp*b3*f0)*(a0*a1 - f0**2) - 2.0_dp*a3*f0*(a0*b1 - b0*f0)
      del3 = (a2*b3 - 3.0_dp*b2*f0)*(a0*a1 - f0**2) + 2.0_dp*f0**2*(3.0_dp*a0*b1 - 2.0_dp*a0*b3 - 3.0_dp*b0*f0)

      !! Find the eta values
      eta0 = del0/delta
      eta1 = del1/delta
      eta2 = del2/delta
      eta3 = del3/delta

      !! Find the p values
      p1p(k) = pi*(1.0_dp + 2.0_dp*R1 + 5.0_dp/4.0_dp*Q1)
      p1m(k) = pi*(1.0_dp - 2.0_dp*R1 + 5.0_dp/4.0_dp*Q1)
      p2p(k) = pi*(1.0_dp + 2.0_dp*R2 + 5.0_dp/4.0_dp*Q2)
      p2m(k) = pi*(1.0_dp - 2.0_dp*R2 + 5.0_dp/4.0_dp*Q2)

      !! Find the q values
      q1p(k) = pi*(-1.0_dp/4.0_dp + 5.0_dp/4.0_dp*Q1 + S1)
      q1m(k) = pi*(-1.0_dp/4.0_dp + 5.0_dp/4.0_dp*Q1 - S1)
      q2p(k) = pi*(-1.0_dp/4.0_dp + 5.0_dp/4.0_dp*Q2 + S2)
      q2m(k) = pi*(-1.0_dp/4.0_dp + 5.0_dp/4.0_dp*Q2 - S2)

      !! Find the Z values
      z1p(k) = pi*(eta0 + 2.0_dp*eta1 + 5.0_dp/4.0_dp*eta2)
      z1m(k) = pi*(eta0 - 2.0_dp*eta1 + 5.0_dp/4.0_dp*eta2)
      z2p(k) = pi*(-1.0_dp/4.0_dp*eta0 + 5.0_dp/4.0_dp*eta2 + 2.0_dp*eta3)
      z2m(k) = pi*(-1.0_dp/4.0_dp*eta0 + 5.0_dp/4.0_dp*eta2 - 2.0_dp*eta3)

      !! Calculate the 4x4 flux matrix elements for this layer
      f(1,1,k) = p1m(k)*e1 ; f(1,2,k) = p1p(k)/e1 ; f(1,3,k) = p2m(k)*e2 ; f(1,4,k) = p2p(k)/e2 
      f(2,1,k) = q1m(k)*e1 ; f(2,2,k) = q1p(k)/e1 ; f(2,3,k) = q2m(k)*e2 ; f(2,4,k) = q2p(k)/e2
      f(3,1,k) = p1p(k)*e1 ; f(3,2,k) = p1m(k)/e1 ; f(3,3,k) = p2p(k)*e2 ; f(3,4,k) = p2m(k)/e2
      f(4,1,k) = q1p(k)*e1 ; f(4,2,k) = q1m(k)/e1 ; f(4,3,k) = q2p(k)*e2 ; f(4,4,k) = q2m(k)/e2

      !! Calculate z up and down level vector following PICASO
      z_u(1,k) = z1m(k) * dir(k+1) ; z_u(2,k) = z2m(k) * dir(k+1) ; z_u(3,k) = z1p(k) * dir(k+1) ; z_u(4,k) = z2p(k) * dir(k+1)
      z_d(1,k) = z1m(k) * dir(k) ; z_d(2,k) = z2m(k) * dir(k) ; z_d(3,k) = z1p(k) * dir(k) ; z_d(4,k) = z2p(k) * dir(k)

    end do

    !! Now add these values to the diagonal matrix representation
    !! Rooney/PICASO uses a banded matrix structure (which seems like the best idea)
    !! PICASO helpfully figured out the banded matrix structure coefficents, so we follow their indexing directly

    !! Upper boundary conditions
    Mb(6,1) = p1m(1)
    Mb(6,2) = q1p(1)
    Mb(5,2) = p1p(1)
    Mb(5,3) = q2m(1)
    Mb(4,3) = p2m(1)
    Mb(4,4) = q2p(1)
    Mb(3,4) = p2p(1)
    Mb(7,1) = q1m(1)

    b_top = 0.0_dp

    BB(1) = b_top - z_d(1,1)
    BB(2) = -b_top/4.0_dp - z_d(2,1)

    !! Lower boundary conditions
    Mb(6,4*nlay-1) = f(3,3,nlay) - w_surf*f(1,3,nlay)
    Mb(6,4*nlay) = f(4,4,nlay) - w_surf*f(2,4,nlay)
    Mb(5,4*nlay) = f(3,4,nlay) - w_surf*f(1,4,nlay)
    Mb(7,4*nlay-2) = f(3,2,nlay) - w_surf*f(1,2,nlay)
    Mb(7,4*nlay-1) = f(4,3,nlay) - w_surf*f(2,3,nlay)
    Mb(8,4*nlay-3) = f(3,1,nlay) - w_surf*f(1,1,nlay)
    Mb(8,4*nlay-2) = f(4,2,nlay) - w_surf*f(2,2,nlay)
    Mb(9,4*nlay-3) = f(4,1,nlay) - w_surf*f(2,1,nlay)

    b_surf = 0.0_dp

    BB(4*nlay-1) = b_surf - z_u(3,nlay) + w_surf*z_u(1,nlay)
    BB(4*nlay) = b_surf - z_u(4,nlay) + w_surf*z_u(2,nlay)

    !! fill remaining rows of matrix

    n = 0
    do i = 1, 4*nlay-4, 4
      Mb(5,2:-4:4) = f02(:-1)
      Mb(5,3:-4:4) = f13(:-1) 
      Mb(4,3:-4:4) = f03(:-1)
      Mb(6,1:-4:4) = f01(:-1)
      Mb(6,2:-4:4) = f12(:-1)
      Mb(6,3:-4:4) = f23(:-1)
      Mb(7,0:-4:4) = f00(:-1)
      Mb(7,1:-4:4) = f11(:-1)
      Mb(7,2:-4:4) = f22(:-1)
      Mb(7,3:-4:4) = f33(:-1)
      Mb(8,0:-4:4) = f10(:-1)
      Mb(8,1:-4:4) = f21(:-1)
      Mb(8,2:-4:4) = f32(:-1)
      Mb(9,0:-4:4) = f20(:-1)
      Mb(9,1:-4:4) = f31(:-1)
      Mb(10,0:-4:4) = f30(:-1)

      B(2:-4:4) = z1mn_down(1:) - z1mn_up(:-1)
      B(3:-4:4) = z2mn_down(1:) - z2mn_up(:-1)
    end do

    do i = 1, lm1
      Mb(5,4::4) = -p1pl(1:)
      Mb(5,5::4) = -q1mn(1:)
      Mb(4,4::4) = -q1mn(1:)
      Mb(4,5::4) = -p1mn(1:)
      Mb(4,6::4) = -q2pl(1:)
        
      Mb(3,4::4) = -p1mn(1:)
      Mb(3,5::4) = -q1pl(1:)
      Mb(3,6::4) = -p2pl(1:)
      Mb(3,7::4) = -q2mn(1:)
        
      Mb(2,5::4) = -p1pl(1:)
      Mb(2,6::4) = -q2mn(1:)
      Mb(2,7::4) = -p2mn(1:)
        
      Mb(1,6::4) = -p2mn(1:)
      Mb(1,7::4) = -q2pl(1:)

      Mb(0,7::4) = -p2pl(1:)
      Mb(6,4::4) = -q1pl(1:)
        
      B(4::4) = z1pl_down(1:) - z1pl_up(:-1)
      B(5::4) = z2pl_down(1:) - z2pl_up(:-1)\
    end do

    !! flux at bottom of atmosphere
    F_bot(4*nlay-3) = f(3,1,nlay)
    F_bot(4*nlay-2) = f(3,2,nlay)
    F_bot(4*nlay-1) = f(3,3,nlay)
    F_bot(4*nlay) = f(3,4,nlay)
    G_bot = z_u(3,nlay)


  end subroutine sw_sph

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

end module ts_spharm_mod
