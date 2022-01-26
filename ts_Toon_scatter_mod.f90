!!!
! Elspeth KH Lee - May 2021
! Two-stream method following the Toon method (Toon et al. 1989)
!
! Based on the CHIMERA code by Mike Line, but cleaned up slightly
! Pros: Fast, accurate at high optical depths, well used and familiar method
! Cons:
! IN DEVELOPMENT - USE WITH UTMOST CAUTION, YOU WILL GET WRONG ANSWERS
!!!

module ts_Toon_scatter_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  real(dp), parameter :: ubari = 0.5_dp

  !! Gauss quadrature variables, cosine angle values (uarr] and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! single angle diffusion factor approximation - typically 1/1.66
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Lacis & Oinas (1991) 3 point numerical values - Does not work somehow, e-mail me if you know why ;)
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

  public :: ts_Toon_scatter
  private :: lw_grey_updown_Toon, sw_grey_updown_Toon, linear_log_interp, bezier_interp

contains

  subroutine ts_Toon_scatter(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB, sw_a_surf, lw_a_surf
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

    !! Output variables
    real(dp),  intent(out) :: olr, asr
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
      call sw_grey_updown_Toon(nlay, nlev, Finc, tau_Ve(:), mu_z, sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_up(:) = 0.0_dp
      sw_down(:) = 0.0_dp
    end if

    ! do i = 1, nlev
    !   print*,i, pe(i)/1e5_dp, sw_up(i), sw_down(i)
    ! end do
    ! stop

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    call lw_grey_updown_Toon(nlay, nlev, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    ! Uncomment if used CHIMERA lower boundary condition
    net_F(nlev) = be_int * pi

    !! Output olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_Toon_scatter

  subroutine lw_grey_updown_Toon(nlay, nlev, be, be_int, tau_IR1, w01, gin, rsurf, lw_up, lw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: rsurf, be_int
    real(dp), dimension(nlev), intent(in) :: tau_IR1, be
    real(dp), dimension(nlay), intent(in) :: w01, gin

    !! Output
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down
    real(dp), dimension(nlev) :: lw_up_raw, lw_down_raw

    !! Work variables
    integer :: k, i, n, m
    integer :: l, lm2, lm1
    real(dp) :: Bsurf, Btop, bottom, tautop
    real(dp), dimension(nlev) :: tau_IR
    real(dp), dimension(nlay) :: dtau1, dtau
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: B0, B1
    real(dp), dimension(nlay) :: lam, gam, alp, term
    real(dp), dimension(nlay) :: Am, Ap, Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xkk
    real(dp), dimension(nlay) :: xk1, xk2

    real(dp), dimension(nlay) :: alphax, g, h, xj, xk
    real(dp), dimension(nlay) :: alpha1, alpha2, sigma1, sigma2
    real(dp), dimension(nlay) :: em1, obj, epp, em4_mp, obj2, epp2, em2, em3

    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l -1

    do k = 1, nlay
      dtau1(k) = max(tau_IR1(k+1) - tau_IR1(k),1e-5_dp)
    end do

    ! Delta eddington scaling
    w0(:) = (1.0_dp - gin(:)**2)*w01(:)/(1.0_dp-w01(:)*gin(:)**2)
    dtau(:) = (1.0_dp-w01(:)*gin(:)**2)*dtau1(:)
    hg(:) = gin(:)/(1.0_dp + gin(:))

    tau_IR(1) = 0.0_dp
    do k = 1, nlay
      tau_IR(k+1) = tau_IR(k) + dtau(k)
    end do

    alp(:) = sqrt((1.0_dp - w0(:))/1.0_dp - w0(:)*hg(:))
    lam(:) = alp(:)*(1.0_dp-w0(:)*hg(:))/ubari
    gam(:) = (1.0_dp-alp(:))/(1.0_dp+alp(:))
    term(:) = ubari/(1.0_dp-w0(:)*hg(:))

    do k = 1, nlay
      if (dtau(k) < 3e-6_dp) then
        ! For low optical depths use the isothermal approimation
        B1(k) = 0.0_dp
        B0(k) = 0.5_dp*(be(k+1) + be(k))
      else
        B1(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
        B0(k) = be(k)
      endif
    end do

    !Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1(:) = B0(:) + B1(:)*term(:)
    Cmm1(:) = B0(:) - B1(:)*term(:)
    !Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp(:) = B0(:) + B1(:)*dtau(:)+B1(:)*term(:)
    Cm(:) = B0(:) + B1(:)*dtau(:)-B1(:)*term(:)

    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp-exp(-tautop/ubari))*be(1)
    Bsurf = be(nlev)
    bottom = Bsurf + B1(nlev)*ubari

    !Solve for the coefficients of system of equations using boundary conditions
    !Exponential terms:
    exptrm(:) = min(lam(:)*dtau(:),35.0_dp)
    Ep(:) = exp(exptrm(:))
    Em(:) = 1.0_dp/Ep(:)

    E1(:) = Ep(:) + gam(:)*Em(:)
    E2(:) = Ep(:) - gam(:)*Em(:)
    E3(:) = gam(:)*Ep(:) + Em(:)
    E4(:) = gam(:)*Ep(:) - Em(:)

    Af(1) = 0.0_dp
    Bf(1) = gam(1) + 1.0_dp
    Cf(1) = gam(1) - 1.0_dp
    Df(1) = btop - Cmm1(1)

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Bf(i) = (E2(n)+E4(n))*(gam(n+1)-1.0_dp)
      Cf(i) = 2.0_dp*(1.0_dp-gam(n+1)**2)
      Df(i) = (gam(n+1)-1.0_dp)*(Cpm1(n+1) - Cp(n)) + (1.0_dp-gam(n+1))*(Cm(n)-Cmm1(n+1))
    end do

    n = 0
    do i = 3, lm1, 2
      n = n + 1
      Af(i) = 2.0_dp*(1.0_dp-gam(n)**2)
      Bf(i) = (E1(n)-E3(n))*(1.0_dp + gam(n+1))
      Cf(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Df(i) = E3(n)*(Cpm1(n+1) - Cp(n)) + E1(n)*(Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - rsurf*E3(nlay)
    Bf(l) = E2(nlay) - rsurf*E4(nlay)
    Cf(l) = 0.0_dp
    Df(l) = bsurf - Cp(nlay) + rsurf*Cm(nlay)

    ! do k = 1, l
    !   print*, k, Af(k), Bf(k), Cf(k), Df(k)
    ! end do
    ! stop

    call dtridgl(l, Af, Bf, Cf, Df, xkk)

    do n = 1, nlay
      xk1(n) = xkk(2*n-1)+xkk(2*n)
      xk2(n) = xkk(2*n-1)-xkk(2*n)
      if (xk2(n) == 0.0_dp) then
        cycle
      end if
      if (abs(xk2(n)/xkk(2*n-1)) < 1e-30_dp) then
        xk2(n) = 0.0_dp
      end if
      !print*, xk1(n), xk2(n)
    end do

    !this would be the raw two stream solution Fluxes, but we
    !won't use these. Will use source function technique
    lw_up_raw(:) = pi*(xk1(:)*Ep(:) + gam(:)*xk2(:)*Em(:) + Cpm1(:))
    lw_down_raw(:) = pi*(xk1(:)*Ep(:)*gam(:) + xk2(:)*Em(:) + Cmm1(:))

    ! do k = 1, nlev
    !   print*, k, lw_up_raw(k), lw_down_raw(k)
    ! end do
    ! stop

    do k = 1, nlay
      if (w0(k) <= 0.01_dp) then
        g(k) = 0.0_dp
        h(k) = 0.0_dp
        xj(k) = 0.0_dp
        xk(k) = 0.0_dp
        alpha1(k)=twopi*B0(k)
        alpha2(k)=twopi*B1(k)
        sigma1(k)=alpha1(k)
        sigma2(k)=alpha2(k)
      else
        alphax(k)=sqrt((1.0_dp-w0(k))/(1.0_dp-w0(k)*hg(k)))
        g(k)=twopi*w0(k)*xk1(k)*(1.0_dp+hg(k)*alphax(k))/(1.0_dp+alphax(k))
        h(k)=twopi*w0(k)*xk2(k)*(1.0_dp-hg(k)*alphax(k))/(1.0_dp+alphax(k))
        xj(k)=twopi*w0(k)*xk1(k)*(1.0_dp-hg(k)*alphax(k))/(1.0_dp+alphax(k))
        xk(k)=twopi*w0(k)*xk2(k)*(1.0_dp+hg(k)*alphax(k))/(1.0_dp+alphax(k))
        alpha1(k)=twopi*(B0(k)+B1(k)*(ubari*w0(k)*hg(k)/(1.0_dp-w0(k)*hg(k))))
        alpha2(k)=twopi*B1(k)
        sigma1(k)=twopi*(B0(k)-B1(k)*(ubari*w0(k)*hg(k)/(1.0_dp-w0(k)*hg(k))))
        sigma2(k)=alpha2(k)
      end if
    end do

    obj(:) = min(lam(:)*dtau(:),35.0_dp)
    em1(:) = exp(-obj(:))
    epp(:) = exp(obj(:))
    obj2(:) = min(0.5_dp*lam(:)*dtau(:),35.0_dp)
    epp2(:) = exp(obj2(:))

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    do m = 1, nmu

      em2(:) = exp(-dtau(:)/uarr(m))
      em3(:) = em1(:)*em2(:)

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = twopi*(1.0_dp - exp(-tautop/uarr(m)))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*em2(k) + &
        & (xj(k)/(lam(k)*uarr(m)+1.0_dp))*(epp(k)-em2(k)) + &
        & (xk(k)/(lam(k)*uarr(m)-1.0_dp))*(em2(k)-em(k))+sigma1(k)*(1.0_dp-em2(k)) + &
        & sigma2(k)*(uarr(m)*em2(k)+dtau(k)-uarr(m))
      end do

      lw_up_g(nlev) = twopi*(Bsurf+B1(nlev)*uarr(m))
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*em2(k) + &
        & (g(k)/(lam(k)*uarr(m)-1.0_dp))*(epp(k)*em2(k)-1.0_dp) + &
        & (h(k)/(lam(k)*uarr(m)+1.0_dp))*(1.0_dp-em3(k))+alpha1(k)*(1.0_dp-em2(k)) + &
        & alpha2(k)*(uarr(m)-(dtau(k)+uarr(m))*em2(k))
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_grey_updown_Toon

  subroutine sw_grey_updown_Toon(nlay, nlev, Finc, tau_V1, mu_z, w01, gin, rsurf, sw_down, sw_up)
    implicit none

    !! Input
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, rsurf
    real(dp), dimension(nlev), intent(in) :: tau_V1
    real(dp), dimension(nlay), intent(in) :: w01, gin

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down

    !! Work variables
    integer :: k, i, n
    integer :: l, lm2, lm1
    real(dp) :: bsurf, btop
    real(dp), dimension(nlev) :: direct, tau_V
    real(dp), dimension(nlay) :: dtau1, dtau
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: g1, g2, g3, g4
    real(dp), dimension(nlay) :: lam, gam, alp, denom
    real(dp), dimension(nlay) :: Am, Ap, Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xk
    real(dp), dimension(nlay) :: xk1, xk2

    !! Optimisation Variables
    real(dp), dimension(nlay) :: opt1

    ! Surface 'emission' boundary fluxes (0 for shortwave)
    bsurf = 0.0_dp
    btop = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w01(:) <= 1.0e-12_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_V1(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - rsurf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l -1

    do k = 1, nlay
      dtau1(k) = max(tau_V1(k+1) - tau_V1(k),1e-5_dp)
    end do

    ! Delta eddington scaling
    w0(:) = (1.0_dp - gin(:)**2)*w01(:)/(1.0_dp-w01(:)*gin(:)**2)
    dtau(:) = (1.0_dp-w01(:)*gin(:)**2)*dtau1(:)
    hg(:) = gin(:)/(1.0_dp + gin(:))

    tau_V(1) = 0.0_dp
    do k = 1, nlay
      tau_V(k+1) = tau_V(k) + dtau(k)
    end do

    direct(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

    g1(:) = 0.86602540378_dp * (2.0_dp-w0(:)*(1.0_dp+hg(:)))
    g2(:) = (1.7320508075688772_dp*w0(:)/2.0_dp) * (1.0_dp-hg(:))
    where (g2(:) == 0.0_dp)
      g2(:) = 1.0e-10_dp
    end where
    g3(:) = (1.0_dp - 1.7320508075688772_dp*hg(:)*mu_z)/2.0_dp
    g4(:) = 1.0_dp - g3(:)

    lam(:) = sqrt(g1(:)**2 - g2(:)**2)
    gam(:) = (g1(:) - lam(:))/g2(:)
    alp(:) = sqrt((1.0_dp - w0(:))/1.0_dp - w0(:)*hg(:))

    denom(:) = lam(:)**2 - 1.0_dp/(mu_z**2)
    where (denom(:) == 0.0_dp)
      denom(:) = 1.0e-10_dp
    end where
    Am(:) = Finc * w0(:) * (g4(:) * (g1(:) + 1.0_dp/mu_z) + g2(:)*g3(:))/denom(:)
    Ap(:) = Finc * w0(:) * (g3(:) * (g1(:) - 1.0_dp/mu_z) + g2(:)*g4(:))/denom(:)

    ! Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    opt1(:) = exp(-tau_V(1:nlay)/mu_z)
    Cpm1(:) = Ap(:) * opt1(:)
    Cmm1(:) = Am(:) * opt1(:)
    ! Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    opt1(:) = exp(-tau_V(2:nlev)/mu_z)
    Cp(:) = Ap(:) * opt1(:)
    Cm(:) = Am(:) * opt1(:)

    ! Solve for the coefficients of system of equations using boundary conditions
    ! Exponential terms:
    exptrm(:) = min(lam(:)*dtau(:),35.0_dp)
    Ep(:) = exp(exptrm(:))
    Em(:) = 1.0_dp/Ep(:)

    E1(:) = Ep(:) + gam(:)*Em(:)
    E2(:) = Ep(:) - gam(:)*Em(:)
    E3(:) = gam(:)*Ep(:) + Em(:)
    E4(:) = gam(:)*Ep(:) - Em(:)

    Af(1) = 0.0_dp
    Bf(1) = gam(1) + 1.0_dp
    Cf(1) = gam(1) - 1.0_dp
    Df(1) = btop - Cmm1(1)

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Bf(i) = (E2(n)+E4(n))*(gam(n+1)-1.0_dp)
      Cf(i) = 2.0_dp*(1.0_dp-gam(n+1)**2)
      Df(i) = (gam(n+1)-1.0_dp)*(Cpm1(n+1) - Cp(n)) + (1.0_dp-gam(n+1))*(Cm(n)-Cmm1(n+1))
    end do

    n = 0
    do i = 3, lm1, 2
      n = n + 1
      Af(i) = 2.0_dp*(1.0_dp-gam(n)**2)
      Bf(i) = (E1(n)-E3(n))*(1.0_dp + gam(n+1))
      Cf(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Df(i) = E3(n)*(Cpm1(n+1) - Cp(n)) + E1(n)*(Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - rsurf*E3(nlay)
    Bf(l) = E2(nlay) - rsurf*E4(nlay)
    Cf(l) = 0.0_dp
    Df(l) = bsurf - Cp(nlay) + rsurf*Cm(nlay)

    call dtridgl(l, Af, Bf, Cf, Df, xk)

    do n = 1, nlay
      xk1(n) = xk(2*n-1)+xk(2*n)
      xk2(n) = xk(2*n-1)-xk(2*n)
      if (xk2(n) == 0.0_dp) then
        cycle
      end if
      if (abs(xk2(n)/xk(2*n-1)) < 1e-30_dp) then
        xk2(n) = 0.0_dp
      end if
    end do

    sw_up(1:nlay) = xk1(:)+gam(:)*xk2(:)+Cpm1(:)
    sw_down(1:nlay) = xk1(:)*gam(:)+xk2(:)+Cmm1(:)

    sw_up(nlev) = xk1(nlay)*Ep(nlay)+gam(nlay)*xk2(nlay)*Em(nlay)+Cp(nlay)
    sw_down(nlev) = xk1(nlay)*Ep(nlay)*gam(nlay)+xk2(nlay)*Em(nlay)+Cm(nlay)

    sw_down(:) = sw_down(:) + direct(:)
    sw_up(:) = sw_up(:)

  end subroutine sw_grey_updown_Toon

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
      !wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      !wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      !if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      !end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      !wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      !wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      !if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      !end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_Toon_scatter_mod
