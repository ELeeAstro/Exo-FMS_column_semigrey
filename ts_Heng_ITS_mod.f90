
!!! Work in progress - do not use!

module ts_Heng_ITS_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !real(dp), parameter :: e2 = 0.5_dp
  real(dp), parameter :: e2 = 2.0_dp/3.0_dp
  !real(dp), parameter :: e2 = 2.0_dp/sqrt(3.0_dp)


  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! single angle diffusion factor approximation - typically 1/1.66
  !integer, parameter :: nmu = 1
  !real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  !real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Lacis & Oinas (1991) 3 point numerical values - Does not work somehow, e-mail me if you know why :)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1_dp, 0.5_dp, 1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = (/0.0433_dp, 0.5742_dp, 0.3825_dp/)

  !! Legendre quadrature for 4 nodes
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.06943184_dp, 0.33000948_dp, 0.66999052_dp, 0.93056816_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.17392742_dp, 0.32607258_dp, 0.32607258_dp, 0.17392742_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! 5 point EGP quadrature values
  ! integer, parameter :: nmu = 5
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   &(/0.0985350858_dp, 0.3045357266_dp, 0.5620251898_dp, 0.8019865821_dp, 0.9601901429_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = &
  !   & (/0.0157479145_dp, 0.0739088701_dp, 0.1463869871_dp, 0.1671746381_dp, 0.0967815902_dp/)

  public :: ts_Heng_ITS
  private :: sw_grey_updown, linear_log_interp, bezier_interp

contains

  subroutine ts_Heng_ITS(Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, sw_a_surf, net_F, olr, asr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, Tint, AB, Ts, sw_a_surf
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g

    !! Output variables
    real(dp), intent(out) :: olr, asr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: lpe
    real(dp), dimension(nlay) :: lTl, lpl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
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
    if (mu_z(nlev) > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_updown(nlay, nlev, Finc, tau_Ve(:), mu_z(:), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature

    !call lw_grey_updown(nlay, nlev, be, be_int, tau_IRe(:), lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output the olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

    stop

  end subroutine ts_Heng_ITS

  subroutine sw_grey_updown(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve, mu_z
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables
    real(dp), dimension(nlay) ::  w, g, f, E
    real(dp), dimension(nlev) :: Fdir
    real(dp), dimension(nlay) :: dtau, T, ccp, ccm
    real(dp), dimension(nlay) :: xi, eps, phi
    real(dp), dimension(nlay) :: Gp, Gm, mu_zm
    real(dp), dimension(nlay) :: Gt1, Gt1_crit, Gt2, Gt3

    !! Matrix variables (Follows Toon variables)
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xk
    real(dp), dimension(nlay) :: xk1, xk2
    real(dp) :: btop, bsurf
    integer :: l, lm2, lm1, n


    w(:) = w_in(:)
    g(:) = g_in(:)

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l - 1

    Fdir(:) = Finc * mu_z(:) * exp(-tau_Ve(:)/mu_z(:))

    E(:) = 1.225_dp - 0.1582_dp*g(:) - 0.1777_dp*w(:) - 0.07465_dp*g(:)**2 + &
      & 0.2351_dp*w(:)*g(:) - 0.05582_dp*w(:)**2

    dtau(:) = tau_Ve(2:) - tau_Ve(1:)

    ccp(:) = 0.5_dp * (1.0_dp + sqrt((E(:)-w(:))/(E(:)*(1.0_dp-w(:)*g(:)))))
    ccm(:) = 0.5_dp * (1.0_dp - sqrt((E(:)-w(:))/(E(:)*(1.0_dp-w(:)*g(:)))))

    T(:) = exp(-2.0_dp * sqrt(E(:)*(E(:)-w(:))*(1.0_dp-w(:)*g(:))) * dtau(:))

    xi(:) = ccm(:)**2*T(:)-ccp(:)**2
    eps(:) = ccp(:)*ccm(:)*(1.0_dp - T(:)**2)
    phi(:) = (ccm(:)**2 - ccp(:)**2)*T(:)

    mu_zm(:) = -(mu_z(1:nlay) + mu_z(2:nlev))/2.0_dp ! Zenith angle at midpoints

    Gt1_crit(:) = 4.0_dp*E(:)*mu_zm(:)**2*(E(:)-w(:))*(1.0_dp-w(:)*g(:))-1.0_dp
    where (Gt1_crit(:) == 0.0_dp)
      Gt1_crit(:) = 2.0_dp*mu_zm(:)**2
    end where

    Gt1(:) = (w(:)*(2.0_dp*E(:)*(1.0_dp-w(:)*g(:)) + g(:)/e2))/(Gt1_crit(:))
    Gt2(:) = 1.0_dp/(2.0_dp*E(:)*(1.0_dp-w(:)*g(:)))
    Gt3(:) = (w(:)*g(:))/(2.0_dp*e2*E(:)*(1.0_dp-w(:)*g(:)))

    Gp(:) = 0.5_dp * (Gt1(:)*(mu_zm(:) + Gt2(:)) + Gt3(:))
    Gm(:) = 0.5_dp * (Gt1(:)*(mu_zm(:) - Gt2(:)) - Gt3(:))

    Ap(:) = (1.0_dp/xi(:)) * (phi(:)*Gp(:)*Fdir(2:nlev) - (eps(:)*Gm(:) + xi(:)*Gp(:))* Fdir(1:nlay))
    Am(:) = (1.0_dp/xi(:)) * (phi(:)*Gm(:)*Fdir(1:nlay) - (eps(:)*Gp(:) + xi(:)*Gm(:))* Fdir(2:nlev))

    ! Surface 'emission' boundary fluxes (0 for shortwave)
    bsurf = 0.0_dp
    btop = 0.0_dp

    Af(1) = 0.0_dp
    Bf(1) =
    Cf(1) =
    Df(1) =

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Bf(i) =
      Cf(i) =
      Df(i) =
    end do



  end subroutine sw_grey_updown

  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    integer, parameter :: nmax = 301
    real(dp) :: x, xkb
    real(dp), dimension(nmax) :: as, ds

    as(l) = af(l)/bf(l)
    ds(l) = df(l)/bf(l)

    do i = l-1, 1, -1
      x = 1.0_dp/(bf(i) - cf(i)*as(i+1))
      as(i) = af(i)*x
      ds(i) = (df(i) - cf(i)*ds(i+1))*x
    end do

    xk(1) = ds(1)

    do i = 2, l
      xk(i) = ds(i)-as(i)*xk(i-1)
    end do

  end subroutine dtridgl

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
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

    real(dp) :: xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

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
      if (w < min(wlim,wlim1) .or. w > max(wlim,wlim1)) then
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
      if (w < min(wlim,wlim1) .or. w > max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_Heng_ITS_mod
