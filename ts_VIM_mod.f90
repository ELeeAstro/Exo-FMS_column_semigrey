!!!
! Elspeth KH Lee - Aug 2023 : Initial version
! sw: Adding layer method with scattering
! lw: 
!     Pros: 
!     Cons:
!!!

module ts_VIM_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Legendre quadrature for 1 nodes
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  public :: ts_VIM
  private :: lw_grey_updown, sw_grey_updown_adding, linear_log_interp, bezier_interp

contains

  subroutine ts_VIM(surf, Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr, net_Fs)
    implicit none

    !! Input variables
    logical, intent(in) :: surf, Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB, sw_a_surf, lw_a_surf, Ts
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
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
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      call sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve(:), mu_z, sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
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
    call lw_grey_updown(surf, nlay, nlev, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

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

  end subroutine ts_VIM

  subroutine lw_grey_updown(surf, nlay, nlev, be, be_int, tau_IRe, ww, gg, lw_a_surf, lw_up, lw_down)
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
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: dtau, beta, eps, dtau_a
    real(dp), dimension(nmu, nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nmu, nmu, nlay) :: Sp, Sm

    real(dp) :: phip, phim, cp, cm, bp, bm, dpp, dm, zepp, zemm, zepm, zemp
    real(dp) :: first, second, third

    
    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:) - tau_IRe(1:nlay)

    ! Delta eddington scaling
    w0(:) = (1.0_dp - gg(:)**2)*ww(:)/(1.0_dp-ww(:)*gg(:)**2)
    dtau(:) = (1.0_dp-ww(:)*gg(:)**2)*dtau(:)
    hg(:) = gg(:)/(1.0_dp + gg(:))

    !! Log B with tau function
    do k = 1, nlay
      if (dtau(k) < 1.0e-7_dp) then
        ! Too low optical depth for numerical stability, Bln = 0
        beta(k) = 0.0_dp
      else
        ! Log B with tau value
        beta(k) = -log(be(k+1)/be(k))/dtau(k)
      end if
    end do

    !! modified co-albedo epsilon
    eps(:) = sqrt((1.0_dp - w0(:))*(1.0_dp - hg(:)*w0(:)))
    !eps(:) = (1.0_dp - w0(:))

    !! Absorption/modified optical depth for transmission function
    dtau_a(:) = eps(:)*dtau(:)

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do i = 1, nmu

      !! Begin two-stream loops
      !! Perform downward loop first - also calculate efficency variables
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(i,1) = 0.0_dp
      do k = 1, nlay

        !! Efficency variables

        !! Downward AA sweep
        lw_down_g(i,k+1) = lw_down_g(i,k)*exp(-dtau_a(k)/uarr(i)) + &
          & (1.0_dp - w0(k))/(uarr(i)*beta(k) - 1.0_dp + w0(k)) * (be(k)*exp(-dtau_a(k)/uarr(i)) - be(k+1))
          !print*, k, m, lw_down_g(k+1), alkap(k), (be(k) - alkap(k)*uarr(m))*exp(-dtau_a(k)/uarr(m))
      end do

      !! Perform upward loop
      if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        lw_up_g(i,nlev) = lw_down_g(i,nlev)*lw_a_surf + be_int
      else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(i,nlev) = lw_down_g(i,nlev) + be_int
      end if
      do k = nlay, 1, -1
        !! Upward AA sweep
        lw_up_g(i,k) = lw_up_g(i,k+1)*exp(-dtau_a(k)/uarr(i)) + &
          & (1.0_dp - w0(k))/(uarr(i)*beta(k) + 1.0_dp - w0(k)) * (be(k) - be(k+1)*exp(-dtau_a(k)/uarr(i)))
      end do

      ! !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      ! lw_down(:) = lw_down(:) + lw_down_g(m,:) * wuarr(m)
      ! lw_up(:) = lw_up(:) + lw_up_g(m,:) * wuarr(m)

    end do

    !! The flux is the integrated intensity * 2pi
    ! lw_down(:) = twopi * lw_down(:)
    ! lw_up(:) = twopi * lw_up(:)
    ! return

    !! Find Sp and Sm - it's now best to put mu into the inner loop
    ! Sp and Sm defined at lower level edges, zero upper boundary condition
    !Sp(:,:) = 0.0_dp
    !Sm(:,:) = 0.0_dp
    do k = 1, nlay
      if (w0(k) <= 1.0e-6_dp) then
         Sp(:,:,k) = 0.0_dp
         Sm(:,:,k) = 0.0_dp
        cycle
      end if
      do i = 1, nmu

        !Sp(i,:,k) = 0.0_dp
        !Sm(i,:,k) = 0.0_dp

        do j = 1, nmu

          phip = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(j)
          phim = 1.0_dp + 3.0_dp*hg(k)*-(uarr(i))*uarr(j)

          !! Note, negative sign mistake in Zhang et al. (2017) - zepm and zepp must be positive quantities
          zepp = -(uarr(i)*uarr(j))/(uarr(i)*(1.0_dp - w0(k)) - uarr(j))
          zemp = (-(uarr(i))*uarr(j))/(-(uarr(i))*(1.0_dp - w0(k)) - uarr(j))
          zepm = -(uarr(i)*-(uarr(j)))/(uarr(i)*(1.0_dp - w0(k)) + uarr(j))
          zemm = (uarr(i)*uarr(j))/(-(uarr(i))*(1.0_dp - w0(k)) + uarr(j))

          cp = (1.0_dp - w0(k))/(uarr(j)*beta(k) + 1.0_dp - w0(k))
          cm = (1.0_dp - w0(k))/(-(uarr(j))*beta(k) + 1.0_dp - w0(k))

          first = phip * zepm * (lw_down_g(j,k) - be(k)*cm) * &
            & (1.0_dp - exp(-dtau(k)/zepm))
          second = phim * zepp * (lw_up_g(j,k+1) - be(k+1)*cp) * &
            & (exp(-(1.0_dp - w0(k))*dtau(k)/uarr(j)) - exp(-dtau(k)/uarr(i)))

          Sp(i,j,k) = first + second

          ! print*, k, i, j
          ! print*, first, second
          ! print*, zepm, uarr(j)/(2.0_dp - w0(k))
          ! print*, zepp, uarr(j)/w0(k)

          first = phip * zemp * (lw_up_g(j,k+1) - be(k+1)*cp) * &
            & (1.0_dp - exp(-dtau(k)/zemp))
          second = phim * zemm * (lw_down_g(j,k) - be(k)*cm) * &
            & (exp(-(1.0_dp - w0(k))*dtau(k)/uarr(j)) - exp(-dtau(k)/uarr(i)))

          Sm(i,j,k) = first + second

          ! print*, k, i, j
          ! print*, zepp, zemp, zepm, zemm
          ! print*, zemp, uarr(j)/(2.0_dp - w0(k))
          ! print*, zemm, uarr(j)/w0(k)

          !print*, k,i,j,Sp(i,j,k)

        end do
      end do
    end do

    !stop

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Do final two sweeps including scattering source function
    do i = 1, nmu

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(i,1) = 0.0_dp
      do k = 1, nlay

        dm = (1.0_dp - w0(k))/(-(uarr(i))*beta(k) + 1.0_dp)

        lw_down_g(i,k+1) = lw_down_g(i,k)*exp(-dtau(k)/uarr(i)) + &
          & dm*(be(k+1) - be(k)*exp(-dtau(k)/uarr(i)))


        bp = uarr(i)/(uarr(i)*beta(k) + 1.0_dp)  
        bm = -(uarr(i))/(-(uarr(i))*beta(k) + 1.0_dp)

        third = 0.0_dp
        do j = 1, nmu

          phip = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(j)
          phim = 1.0_dp + 3.0_dp*hg(k)*-(uarr(i))*uarr(j)

          cp = (1.0_dp - w0(k))/(uarr(j)*beta(k) + 1.0_dp - w0(k))
          cm = (1.0_dp - w0(k))/(-(uarr(j))*beta(k) + 1.0_dp - w0(k))

          third = third + &
            & (Sm(i,j,k) - bm*(cp*phim + cm*phip)*(be(k+1) - be(k)*exp(-dtau(k)/uarr(i))))
          !print*, k, m, lw_down_g(k+1), alkap(k), (be(k) - alkap(k)*uarr(m))*exp(-dtau_a(k)/uarr(m))

        end do

        !print*, k,third

        lw_down_g(i,k+1) = lw_down_g(i,k+1) - w0(k)/(real(nmu*2,dp)*(-uarr(i)))*third
      end do

      !! Perform upward loop
      if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        lw_up_g(i,nlev) = lw_down_g(i,nlev)*lw_a_surf + be_int
      else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(i,nlev) = lw_down_g(i,nlev) + be_int
      end if

      do k = nlay, 1, -1

        dpp = (1.0_dp - w0(k))/(uarr(i)*beta(k) + 1.0_dp)

        lw_up_g(i,k) = lw_up_g(i,k+1)*exp(-dtau(k)/uarr(i)) - &
          & dpp*(be(k+1)*exp(-dtau(k)/uarr(i)) - be(k))


        bp = uarr(i)/(uarr(i)*beta(k) + 1.0_dp)  

        third = 0.0_dp
        do j = 1, nmu

          phip = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(j)
          phim = 1.0_dp + 3.0_dp*hg(k)*-(uarr(i))*uarr(j)

          cp = (1.0_dp - w0(k))/(uarr(j)*beta(k) + 1.0_dp - w0(k))
          cm = (1.0_dp - w0(k))/(-(uarr(j))*beta(k) + 1.0_dp - w0(k))

          third = third + &
            & (Sp(i,j,k) - bp*(cp*phip + cm*phim)*(be(k+1)*exp(-dtau(k)/uarr(i)) - be(k)))
          !print*, k, m, lw_down_g(k+1), alkap(k), (be(k) - alkap(k)*uarr(m))*exp(-dtau_a(k)/uarr(m))

        end do

        !print*, k, third

        lw_up_g(i,k) = lw_up_g(i,k) + w0(k)/(real(nmu*2,dp)*(uarr(i)))*third

      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(i,:) * wuarr(i)
      lw_up(:) = lw_up(:) + lw_up_g(i,:) * wuarr(i)

    end do

    !! The flux is the integrated intensity * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)

  end subroutine lw_grey_updown


  subroutine sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
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

    ! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
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
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_z**2)/(1.0_dp - lam(k)**2*mu_z**2)
      alp(k) = 0.75_dp * w_s(k) * mu_z * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_z**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_Ve_s(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_Ve_s(k)/mu_z,99.0_dp)
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
    sw_down(:) = sw_down(:) * mu_z * Finc
    sw_up(:) = sw_up(:) * mu_z * Finc

  end subroutine sw_grey_updown_adding

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

  end subroutine bezier_interp

end module ts_VIM_mod
