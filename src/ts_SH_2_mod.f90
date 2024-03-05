!!!
! Elspeth K.H. Lee - 
! sw: 
! lw: 
! Pros:
! Cons:
!!!

module ts_SH_2_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! single angle diffusion factor approximation - typically 1/1.66
  integer, parameter :: nmu = 1
  real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Use Two-Term HG function for sw or lw
  logical, parameter :: TTHG_sw = .False.
  logical, parameter :: TTHG_lw = .False.

  public :: ts_SH_2
  private :: sw_SH_2, lw_SH_2, linear_log_interp, bezier_interp

contains

  subroutine ts_SH_2(surf, Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
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
      call sw_SH_2(nlay, nlev, Finc, tau_Ve(:), mu_z(nlev), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
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

    call lw_SH_2(surf, nlay, nlev, pe, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Net surface flux (for surface temperature evolution)
    !! We have to define positive as downward (heating) and cooling (upward) in this case
    net_Fs = sw_down(nlev) + lw_down(nlev) - lw_up(nlev)

    ! Uncomment if used CHIMERA lower boundary condition
    net_F(nlev) = be_int * pi

    !! Output the olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_SH_2

  subroutine lw_SH_2(surf, nlay, nlev, pe, be, be_int, tau_IRe, w_in, g_in, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe, pe
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: l, k, i, m
    integer, parameter :: nstr = 2
    real(dp), dimension(nstr, nlay) :: a, bb, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface

    real(dp), dimension(2,nlay) :: eta
    real(dp), dimension(nlay) :: lam, expo, exptrm, q, Q1, Q2
    real(dp), dimension(nlay) :: Q1mn, Q2mn, Q1pl, Q2pl
    real(dp), dimension(nlay) :: zmn_up, zpl_up, zmn_down, zpl_down
    real(dp), dimension(nlev) :: tau_e

    real(dp), parameter :: mu1 = 0.5_dp

    real(dp), dimension(nlay) :: b0
    real(dp), dimension(nlay) :: b1
    real(dp) :: tau_top, b_top

    real(dp), dimension(5,2*nlay) :: Mb
    real(dp), dimension(7,2*nlay) :: Mb_F
    real(dp), dimension(2*nlay) :: B

    integer :: info
    integer, dimension(2*nlay) :: ipiv

    real(dp), dimension(2) :: Pubar1
    real(dp), dimension(nlay) :: X1, X2, alpha, beta, expo_alp, expo_bet, exptrm_alp, exptrm_bet
    real(dp), dimension(nlay) :: expdtau, Aint0, Aint1, Nint0, Nint1
    real(dp), dimension(nlay) :: multi_scat, intgrl_per_layer

    real(dp), dimension(nlev) :: lw_up_g, lw_down_g


    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, g0


    ! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:nlev) - tau_IRe(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling (g**nstream)
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

    g0(:) = g_in(:)

    !! Reform edge optical depths
    tau_e(1) = tau_IRe(1)
    do k = 1, nlay
      tau_e(k+1) = tau_e(k) + dtau(k)
    end do

    !! Linear B with tau function
    where (dtau(:) < 1.0e-6_dp)
      b1(:) = 0.0_dp
      b0(:) = 0.5_dp*(be(2:nlev) + be(1:nlay))
    elsewhere
      b1(:) = (be(2:nlev) - be(1:nlay))/dtau(:) 
      b0(:) = be(1:nlay)
    end where

    tau_top = dtau(1)*pe(1)/(pe(2)-pe(1))
    b_top = pi*(1.0_dp - exp(-tau_top / mu1 )) * be(1)

    b_surface = pi*(be(nlev) + b1(nlay)*mu1)
    surf_reflect = 0.0_dp

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (g0(:) - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:)
      bb(l,:) = 0.0_dp
    end do

    eta(1,:) = 0.0_dp
    eta(2,:) = 0.0_dp

    lam(:) = sqrt(a(1,:)*a(2,:))
    expo(:) = lam(:)*dtau(:)
    where (expo(:) > 35.0_dp) 
      expo(:) = 35.0_dp
    end where
    exptrm(:) = exp(-expo(:))

    q(:) = lam(:)/a(2,:)
    Q1(:) = (0.5_dp + q(:))*2.0_dp*pi
    Q2(:) = (0.5_dp - q(:))*2.0_dp*pi

    Q1mn(:) = Q1(:)*exptrm(:);  Q2mn(:) = Q2(:)*exptrm(:)
    Q1pl(:) = Q1(:)/exptrm(:);  Q2pl(:) = Q2(:)/exptrm(:)

    zmn_down(:) = ((1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp - B1(:)/a(2,:))) *2.0_dp*pi           
    zmn_up(:) = ((1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp - B1(:)/a(2,:) + B1*dtau(:)/2.0_dp)) *2.0_dp*pi 
    zpl_down(:) = ((1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp + B1(:)/a(2,:))) *2.0_dp*pi           
    zpl_up(:) = ((1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp + B1(:)/a(2,:) + B1*dtau(:)/2.0_dp)) *2.0_dp*pi 

    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    Mb(3,1) = Q1(1)
    Mb(2,2) = Q2(1)
    B(1) = b_top - zmn_down(1)

    Mb(4, 2*nlay-1) = Q2mn(nlay) - surf_reflect*Q1mn(nlay)
    Mb(3, 2*nlay) = Q1pl(nlay) - surf_reflect*Q2pl(nlay)
    B(2*nlay) = b_surface - zpl_up(nlay) + surf_reflect*zmn_up(nlay)

    do i = 2, nlay
      Mb(1, 2*i) = -Q2(i)
      Mb(2, 2*i - 1) = -Q1(i)
      Mb(2, 2*i) = -Q1(i)
      Mb(3, 2*i - 1) = -Q2(i)
    end do

    do i = 1, nlay-1
      Mb(3, 2*i) = Q2pl(i)
      Mb(4, 2*i-1) = Q1mn(i)
      Mb(4, 2*i) = Q1pl(i)
      Mb(5, 2*i-1) = Q2mn(i)
    end do

    do i = 1, nlay-1
      B(2*i) = zmn_down(i+1) - zmn_up(i)
      B(2*i + 1) = zpl_down(i+1) - zpl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(3,:) = Mb(1,:)
    Mb_F(4,:) = Mb(2,:)
    Mb_F(5,:) = Mb(3,:)  
    Mb_F(6,:) = Mb(4,:)
    Mb_F(7,:) = Mb(5,:)   

    call dgbsv(2*nlay, 2, 2, 1, Mb_F, 7, ipiv, B, 2*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

    !! Split B into downward and upward components
    do i = 1, nlay
      X1(i) = B(2*i - 1)
      X2(i) = B(2*i)
    end do

    !! Now do the source function technique to get up and down fluxes
    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    do m = 1, nmu

      Pubar1(1) = 1.0_dp 
      Pubar1(2) = -uarr(m)

      alpha(:) = 1.0_dp/uarr(m) + lam(:)
      beta(:) = 1.0_dp/uarr(m) - lam(:)
      expo_alp(:) = min(alpha(:) * dtau(:), 35.0)
      expo_bet(:) = min(beta(:) * dtau(:), 35.0) 
      exptrm_alp(:) = (1.0_dp - exp(-expo_alp(:))) / alpha(:)
      exptrm_bet(:) = (1.0_dp - exp(-expo_bet(:))) / beta(:)

      Aint0(:) = X1(:)  * (w_multi(1,:)-w_multi(2,:)*Pubar1(2)*q(:)) * exptrm_alp(:)
      Aint1(:) = X2(:) * (w_multi(1,:)+w_multi(2,:)*Pubar1(2)*q(:)) * exptrm_bet(:)

      expdtau(:) = exp(-dtau(:)/uarr(m))
      Nint0(:) = w_multi(1,:)*((1.0_dp-w0(:)) * uarr(m) & 
        & / a(1,:) * (b0(:) *(1.0_dp-expdtau(:)) + b1(:)*(uarr(m) - (dtau(:)+uarr(m))*expdtau(:)))) 
      Nint1(:) = w_multi(2,:)*Pubar1(2)*((1.0_dp-w0(:)) * uarr(m) / a(1,:) * ( b1(:)*(1-expdtau(:)) / a(2,:)))

      multi_scat(:) = Aint0(:) + Nint0(:) + Aint1(:) + Nint1(:)

      expo(:) = min(dtau(:) / uarr(m), 35.0)  
      expdtau(:) = exp(-expo(:))

      intgrl_per_layer(:) = (w0(:) *  multi_scat(:) *2*pi &
        & + 2.0_dp*pi*(1.0_dp-w0(:)) * uarr(m) * &
        & (b0(:) * (1.0_dp - expdtau(:)) &
        & + b1(:) * (uarr(m) - (dtau(:) + uarr(m)) * expdtau(:))))

      lw_down_g(1) = 2.0_dp*pi*(1.0_dp - exp(-tau_top/uarr(m)))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = (lw_down_g(k) * exp(-dtau(k)/uarr(m)) + intgrl_per_layer(k) / uarr(m)) 
      end do

      lw_up_g(nlev) = be(nlev) + b1(nlay) * uarr(m)*2*pi
      do k = nlay, 1, -1
        lw_up_g(k) = (lw_up_g(k+1) * exp(-dtau(k)/uarr(m)) + intgrl_per_layer(k) / uarr(m)) 
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_SH_2

  subroutine sw_SH_2(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables and arrays
    integer :: l, k, i
    integer, parameter :: nstr = 2
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nstr, nlay) :: a, bb, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface, b_top

    real(dp), dimension(2,nlay) :: eta
    real(dp), dimension(nlay) :: del, lam, expo, exptrm, q, Q1, Q2
    real(dp), dimension(nlay) :: Q1mn, Q2mn, Q1pl, Q2pl
    real(dp), dimension(nlay) :: zmn, zpl, zmn_up, zpl_up, zmn_down, zpl_down
    real(dp), dimension(nlev) :: expon, tau_e

    real(dp), dimension(5,2*nlay) :: Mb
    real(dp), dimension(7,2*nlay) :: Mb_F
    real(dp), dimension(2*nlay) :: B, X

    real(dp), dimension(2*nlev, 2*nlay) :: F
    real(dp), dimension(2*nlev) :: G
    real(dp), dimension(2*nlay) :: F_bot
    real(dp) :: G_bot

    real(dp), dimension(2*nlev) :: flux_temp
    real(dp) :: flux_bot


    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, g0
    integer :: info
    integer, dimension(2*nlay) :: ipiv

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-6_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_Ve(2:nlev) - tau_Ve(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling (g**nstream)
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

    g0(:) = g_in(:)

    !! Reform edge optical depths
    tau_e(1) = tau_Ve(1)
    do k = 1, nlay
      tau_e(k+1) = tau_e(k) + dtau(k)
    end do

    Pu0(1) = 1.0_dp
    Pu0(2) = -mu_z

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (g0(:) - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:)
      bb(l,:) = ((w0(:) * w_multi(l,:)) * Finc * Pu0(l)) / (4.0_dp*pi)
    end do

    surf_reflect = 0.0_dp

    b_surface = 0.0_dp + surf_reflect*mu_z*Finc*exp(-tau_e(nlev)/mu_z)

    b_top = 0.0_dp

    lam(:) = sqrt(a(1,:)*a(2,:))
    expo(:) = min(lam(:)*dtau(:),35.0_dp)
    exptrm(:) = exp(-expo(:))

    del(:) = 1.0_dp/((1.0_dp / mu_z**2) - lam(:)**2)
    eta(1,:) = (bb(2,:)/mu_z - a(2,:)*bb(1,:)) * del(:)
    eta(2,:) = (bb(1,:)/mu_z - a(1,:)*bb(2,:)) * del(:)

    q(:) = lam(:)/a(2,:)
    Q1(:) = (1.0_dp + 2.0_dp*q(:))*pi
    Q2(:) = (1.0_dp - 2.0_dp*q(:))*pi

    Q1mn(:) = Q1(:)*exptrm(:);  Q2mn(:) = Q2(:)*exptrm(:)
    Q1pl(:) = Q1(:)/exptrm(:);  Q2pl(:) = Q2(:)/exptrm(:)

    zmn(:) = (eta(1,:) - 2.0_dp*eta(2,:))*pi
    zpl(:) = (eta(1,:) + 2.0_dp*eta(2,:))*pi
    expon(:) = exp(-tau_e(:)/mu_z)
    zmn_up(:) = zmn(:) * expon(2:nlev)
    zpl_up(:) = zpl(:) * expon(2:nlev) 
    zmn_down(:) = zmn(:) * expon(1:nlay)
    zpl_down(:) = zpl(:) * expon(1:nlay) 

    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    Mb(3,1) = Q1(1)
    Mb(2,2) = Q2(1)
    B(1) = b_top - zmn_down(1)

    Mb(4, 2*nlay-1) = Q2mn(nlay) - surf_reflect*Q1mn(nlay)
    Mb(3, 2*nlay) = Q1pl(nlay) - surf_reflect*Q2pl(nlay)
    B(2*nlay) = b_surface - zpl_up(nlay) + surf_reflect*zmn_up(nlay)

    do i = 2, nlay
      Mb(1, 2*i) = -Q2(i)
      Mb(2, 2*i - 1) = -Q1(i)
      Mb(2, 2*i) = -Q1(i)
      Mb(3, 2*i - 1) = -Q2(i)
    end do

    do i = 1, nlay-1
      Mb(3, 2*i) = Q2pl(i)
      Mb(4, 2*i-1) = Q1mn(i)
      Mb(4, 2*i) = Q1pl(i)
      Mb(5, 2*i-1) = Q2mn(i)
    end do

    do i = 1, nlay-1
      B(2*i) = zmn_down(i+1) - zmn_up(i)
      B(2*i + 1) = zpl_down(i+1) - zpl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(3,:) = Mb(1,:)
    Mb_F(4,:) = Mb(2,:)
    Mb_F(5,:) = Mb(3,:)  
    Mb_F(6,:) = Mb(4,:)
    Mb_F(7,:) = Mb(5,:)   

    call dgbsv(2*nlay, 2, 2, 1, Mb_F, 7, ipiv, B, 2*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

    ! flux at bottom of atmosphere
    F_bot(2*nlay - 1) = Q2mn(nlay)
    F_bot(2*nlay) = Q1pl(nlay)
    G_bot = zpl_up(nlay)

    F(1, 1) = Q1(1)
    F(1, 2) = Q2(1)
    F(2, 1) = Q2(1)
    F(2, 2) = Q1(1)

    k = 0
    do i = 1, 2*nlay, 2
      F(i + 2, i) = Q1mn(k + 1)
      F(i + 2, i + 1) = Q2pl(k + 1)
      F(i + 3, i) = Q2mn(k + 1)
      F(i + 3, i + 1) = Q1pl(k + 1)
      k = k + 1
    end do

    G(1) = zmn_down(1)
    G(2) = zpl_down(1)

    do i = 1, nlay
      G(i*2 + 1) = zmn_up(i)
    end do

    do i = 2, nlay
      G(i*2) = zpl_up(i)
    end do

    flux_temp(:) = matmul(F(:,:),B(:)) + G(:)
    flux_bot = sum(F_bot(:)*X(:)) + G_bot

    do i = 1, nlev
      sw_down(i) = max(flux_temp(i*2-1),0.0_dp) + mu_z * Finc * expon(i)
      sw_up(i) = max(flux_temp(i*2), 0.0_dp)
    end do

    !! Now add these values to the diagonal matrix representation
    !! Rooney/PICASO uses a banded matrix structure (which seems like the best idea)
    !! PICASO helpfully figured out the banded matrix structure coefficents, so we follow their indexing directly
    !! The main banded Matix (Mb) is 5x(2*nlay), 2 streams, so 5 terms needed and 2*nlay for the 2 Flux+ terms
    

  end subroutine sw_SH_2

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

end module ts_SH_2_mod
