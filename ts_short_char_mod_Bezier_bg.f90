!!!
! Elspeth KH Lee - Jun 2022 : Initial version
!
! sw: Adding layer method with scattering
! lw: Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
!     Uses the method of short characteristics (Olson & Kunasz 1987) with Bezier interpolants (de la Cruz Rodriguez and Piskunov 2013 [CR&P13]).
!     Pros: Very fast, accurate at high optical depths, very stable, extreamly smooth heating profiles
!     Cons: No lw scattering
!!!

module ts_short_char_mod_Bezier_bg
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

  public :: ts_short_char_Bezier_bg
  private :: bg_thermal_updown_Bezier, bg_stellar_updown_adding, linear_interp, linear_log_interp, calc_cff,&
             Bezier_interp, Bezier_interp_yc, BB_integrate

contains

  subroutine ts_short_char_Bezier_bg(Bezier, surf, Planck_table, n_bands, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, tau_bg, &
    & mu_z, F0, Tint, AB, sw_a, sw_g, sw_a_surf, lw_a_surf, net_F, thermal_up, thermal_down, stellar_up, stellar_down, opr, asr, &
    & net_Fs, cff, scff, f_star)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier, surf
    integer, intent(in) :: n_bands, nlay, nlev
    real(dp), intent(in) :: Ts, F0, Tint, AB, sw_a_surf, lw_a_surf
    real(dp), dimension(n_bands), intent(in) :: f_star
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe, mu_z
    real(dp), dimension(n_bands+1,2999), intent(in) :: Planck_table
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(n_bands,nlev), intent(in) :: tau_bg
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g

    !! Output variables
    real(dp), dimension(n_bands), intent(out) :: opr, asr, net_Fs, scff
    real(dp), dimension(n_bands,nlev), intent(out) :: net_F
    real(dp), dimension(n_bands,nlev), intent(out) :: thermal_up, thermal_down, stellar_up, stellar_down
    real(dp), dimension(n_bands,nlay), intent(out) :: cff

    !! Work variables
    integer :: i, b, T_closest
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: lpe
    real(dp), dimension(nlay) :: lTl, lpl
    real(dp), dimension(n_bands,nlev) :: be_bg, thermal_net, stellar_net
    real(dp), dimension(2999) :: T_range
    real(dp), dimension(n_bands,2) :: wn_edges
    logical :: integration

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
      if (n_bands .eq. 1) then
        Finc = (1.0_dp - AB) * F0
        call bg_stellar_updown_adding(nlay, nlev, Finc, tau_bg(4,:), mu_z(:), sw_a, sw_g, sw_a_surf, stellar_down(1,:),&
                                    & stellar_up(1,:))
      else if (n_bands .ge. 4) then
        do b = 1, n_bands
          Finc = (1.0_dp - AB) * F0 * f_star(b)
          call bg_stellar_updown_adding(nlay, nlev, Finc, tau_bg(b,:), mu_z(:), sw_a, sw_g, sw_a_surf, stellar_down(b,:),&
                                      & stellar_up(b,:))
        end do
      end if
    else
      stellar_down(:,:) = 0.0_dp
      stellar_up(:,:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    if (surf .eqv. .True.) then
      be_int = (sb * Ts**4)/pi
    else
      be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    end if

    if (n_bands .eq. 1) then

      !be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
      T_range = Planck_table(1,:)
      do i = 1,nlev
        ! be cannot be zero. lw_grey_updown_linear divides by be(k) to get Am(k)
        T_closest = minloc(abs(T_range-(T_range/T_range)*Te(i)), DIM=1)
        be(i) = max(Planck_table(2,T_closest), 1d-30)/pi
      end do

      call bg_thermal_updown_Bezier(surf, nlay, nlev, be, be_int, lw_a_surf, tau_bg(1,:), thermal_up(1,:), thermal_down(1,:))

      !! Net fluxes at each level
      thermal_net(1,:) = thermal_up(1,:) - thermal_down(1,:)
      stellar_net(1,:) = stellar_up(1,:) - stellar_down(1,:)
      net_F(1,:) = thermal_net(1,:) + thermal_net(1,:)

      !! Net surface flux (for surface temperature evolution)
      !! We have to define positive as downward (heating) and cooling (upward) in this case
      net_Fs = thermal_down(1,nlev) + thermal_down(1,nlev) - thermal_up(1,nlev)

      !! Output the opr
      opr(1) = thermal_up(1,1)

      !! Output asr
      asr(1) = stellar_down(1,1) - stellar_up(1,1)

      call calc_cff(surf, nlay, be, pe, tau_bg(1,:), cff(1,:), scff(1))


    else if (n_bands .ge. 4) then

      integration = .false.
      if (integration .eqv. .true.) then
        if (n_bands .eq. 4) then
          wn_edges(1,1) = 12500.0_dp; wn_edges(1,2) = 1.0_dp     ! Region IR + W1 + W2 (requires subtracting B(W1)+B(W2) from B(IR))
          wn_edges(2,1) = 2900.0_dp;  wn_edges(2,2) = 2200.0_dp  ! Region W1
          wn_edges(3,1) = 1300.0_dp;  wn_edges(3,2) = 500.0_dp   ! Region W2
          wn_edges(4,1) = 40500.0_dp; wn_edges(4,2) = 12500.0_dp ! Region SW

        else if (n_bands .eq. 5) then
          wn_edges(1,1) = 12500.0_dp; wn_edges(1,2) = 1.0_dp     ! Region IR + W1 + W2
          wn_edges(2,1) = 2900.0_dp;  wn_edges(2,2) = 2200.0_dp  ! Region W1
          wn_edges(3,1) = 1300.0_dp;  wn_edges(3,2) = 500.0_dp   ! Region W2
          wn_edges(4,1) = 40500.0_dp; wn_edges(4,2) = 25000.0_dp ! Region UV
          wn_edges(5,1) = 25000.0_dp; wn_edges(5,2) = 12500.0_dp ! Region VIS

        else if (n_bands .eq. 6) then
          wn_edges(1,1) = 12500.0_dp; wn_edges(1,2) = 1.0_dp     ! Region IR + W1 + W2real(dp), dimension(n_bands,nlay), intent(out) :: cff
          wn_edges(2,1) = 2900.0_dp;  wn_edges(2,2) = 2200.0_dp  ! Region W1
          wn_edges(3,1) = 1300.0_dp;  wn_edges(3,2) = 500.0_dp   ! Region W2
          wn_edges(4,1) = 40500.0_dp; wn_edges(4,2) = 25000.0_dp ! Region UV
          wn_edges(5,1) = 25000.0_dp; wn_edges(5,2) = 18750.0_dp ! Region VIS1
          wn_edges(6,1) = 18750.0_dp; wn_edges(6,2) = 12500.0_dp ! Region VIS2

        end if

        do b = 1, n_bands
          call BB_integrate(nlev, Te(:), wn_edges(b,:), be_bg(b,:))
        end do

      else
        T_range = Planck_table(1,:)
        !! Fill
        do i = 1,nlev
          ! be cannot be zero. lw_grey_updown_linear divides by be(k) to get Am(k)
          T_closest = minloc(abs(T_range-(T_range/T_range)*Te(i)), DIM=1) ! returns the index of the closest element to Te(i) in TL
          be_bg(1,i) = max(Planck_table(2,T_closest), 1d-30)
          be_bg(2,i) = max(Planck_table(3,T_closest), 1d-30)
          be_bg(3,i) = max(Planck_table(4,T_closest), 1d-30)
          if (n_bands .eq. 4) then
            be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
          elseif (n_bands .eq. 5) then
            be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
            be_bg(5,i) = max(Planck_table(6,T_closest), 1d-30)
          elseif (n_bands .eq. 6) then
            be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
            be_bg(5,i) = max(Planck_table(6,T_closest), 1d-30)
            be_bg(6,i) = max(Planck_table(7,T_closest), 1d-30)
          end if
        end do
      end if

      be_bg = be_bg/pi

      do b = 1, n_bands

        call bg_thermal_updown_Bezier(surf, nlay, nlev, be_bg(b,:), be_int, lw_a_surf, tau_bg(b,:), thermal_up(b,:),&
                                    & thermal_down(b,:))

        !! Net fluxes at each level
        thermal_net(b,:) = thermal_up(b,:) - thermal_down(b,:)
        stellar_net(b,:) = stellar_up(b,:) - stellar_down(b,:)
        net_F(b,:) = thermal_net(b,:) + thermal_net(b,:)

        !! Net surface flux (for surface temperature evolution)
        !! We have to define positive as downward (heating) and cooling (upward) in this case
        net_Fs(b) = thermal_down(b,nlev) + thermal_down(b,nlev) - thermal_up(b,nlev)

        !! Output the opr
        opr(b) = thermal_up(b,1)

        !! Output asr
        asr(b) = stellar_down(b,1) - stellar_up(b,1)

        call calc_cff(surf, nlay, be_bg(b,:), pe, tau_bg(b,:), cff(b,:), scff(b))
      end do

    end if

  end subroutine ts_short_char_Bezier_bg

  subroutine bg_thermal_updown_Bezier(surf, nlay, nlev, be, be_int, lw_a_surf, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, del, edel, del2, del3, edel_lim
    real(dp), dimension(nlay) :: tau_mid, ltau_mid
    real(dp), dimension(nlev) :: lbe, ltau
    real(dp), dimension(nlay) :: Ak, Bk, Gk, Ck
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:) - tau_IRe(1:nlay)
    tau_mid(:) = (tau_IRe(2:) + tau_IRe(1:nlay))/2.0_dp

    ltau(:) = log10(tau_IRe(:))
    ltau_mid(:) = log10(tau_mid(:))
    lbe(:) = log10(be(:))

    ! Find the source function at Bezier control point and the center of optical depth space of each layer
    call bezier_interp_yc(ltau(1:3), lbe(1:3), 3, ltau_mid(1), Ck(1))
    !print*, "ltau(1:3), lbe(1:3), ltau_mid(1), Ck(1) = ", ltau(1:3), lbe(1:3), ltau_mid(1), Ck(1)
    do k = 2, nlay
      call bezier_interp_yc(ltau(k-1:k+1), lbe(k-1:k+1), 3, ltau_mid(k), Ck(k))
      !print*, "ltau(k-1:k+1), lbe(k-1:k+1), ltau_mid(k), Ck(k) = ", ltau(k-1:k+1), lbe(k-1:k+1), ltau_mid(k), Ck(k)
    end do
    Ck(:) = 10.0_dp**Ck(:)

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)/uarr(m)
      del2(:) = del(:)**2
      del3(:) = del(:)**3
      edel(:) = exp(-del(:))
      edel_lim(:) = &
      & max(1.0_dp - del(:) + del2(:)/2.0_dp - del3(:)/6.0_dp, 0.0_dp)/edel(:)

      !! Prepare loop
      ! de la Cruz Rodriguez and Piskunov 2013 Bezier interpolant parameters
      where (edel_lim(:) > 0.999_dp)
        ! If we are in very low optical depth regime,
        ! then use a Taylor expansion following [CR&P13]
        Ak(:) = del(:)/3.0_dp - del2(:)/12.0_dp + del3(:)/60.0_dp
        Bk(:) = del(:)/3.0_dp - del2(:)/4.0_dp + del3(:)/10.0_dp
        Gk(:) = del(:)/3.0_dp - del2(:)/6.0_dp + del3(:)/20.0_dp
      elsewhere
        ! Use Bezier interpolants
        Ak(:) = (2.0_dp + del2(:) - 2.0_dp*del(:) - 2.0_dp*edel(:))/del2(:)
        Bk(:) = (2.0_dp - (2.0_dp + 2.0_dp*del(:) + del2(:))*edel(:))/del2(:)
        Gk(:) = (2.0_dp*del(:) - 4.0_dp + (2.0_dp*del(:) + 4.0_dp)*edel(:))/del2(:)
      end where

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
         lw_down_g(k+1) = lw_down_g(k)*edel(k) + Ak(k)*be(k+1) + Bk(k)*be(k) + Gk(k)*Ck(k)! TS intensity
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
        lw_up_g(k) = lw_up_g(k+1)*edel(k)  + Ak(k)*be(k) + Bk(k)*be(k+1) + Gk(k)*Ck(k)! TS intensity
        !print*, k, lw_up_g(k), lw_down_g(k), Ak(k), Bk(k), Gk(k)
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do
    !! The flux is the intensity * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)

  end subroutine bg_thermal_updown_Bezier

  subroutine bg_stellar_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
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

  end subroutine bg_stellar_updown_adding

  subroutine calc_cff(surf, nlay, be, pe, tau, cff, scff )
    implicit none
    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay
    real(dp), dimension(nlay+1), intent(in) :: be, pe, tau

    !! Output variables
    real(dp), dimension(nlay), intent(out) :: cff
    real(dp), intent(out) :: scff

    !! Work variables and arrays
    integer :: i
    real(dp) :: D_factor

    ! Contribution function to the planetary outgoing radiation
    D_factor = 1.66_dp ! constant diffusivity factor
    do i = 1, nlay
      !D_factor = 1.5_dp + 0.5_dp/(1_dp + 4_dp*(tau_struc_edges(k)+tau_struc_edges(k+1))/2_dp &
      !         + 10_dp*((tau_struc_edges(k)+tau_struc_edges(k+1))/2_dp)**2) ! diffusivity approximation (Thomas & Stammes 1999, section 11.2.5)
      cff(i) = ( 2.d0*pi*be(i)/D_factor ) * &
      & abs(exp(-D_factor*tau(i+1)) - exp(-D_factor*tau(i)) ) / &
      & abs( log10(pe(i+1))-log10(pe(i)) )
    end do
    !scff = 2_dp*pi*flux_up(nlay1)*exp(-D_factor*tau_struc_edges(nlay1))/D_factor
    scff = be(nlay+1)*exp(-D_factor*tau(nlay+1))

  end subroutine calc_cff

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

  ! Perform linear interpolation
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  subroutine Bezier_interp_yc(xi, yi, ni, x, yc)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: yc

    real(dp) :: dx, dx1, dy, dy1, w, wlim, wlim1

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
    end if

  end subroutine Bezier_interp_yc

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

  subroutine BB_integrate(nlay1, Tint, wn_edges, bint)
    implicit none

    integer, intent(in) :: nlay1
    real(kind=dp), dimension(nlay1), intent(in) :: Tint
    real(kind=dp), dimension(2), intent(in) :: wn_edges
    real(kind=dp), dimension(nlay1), intent(out) :: bint

    integer :: i, w, j, intitera

    real(kind=dp), dimension(2) :: iB

    real(kind=dp) :: c1, x, x2, x3, itera, summ, n2, dn
    real(kind=dp) :: h, c, k, c2

    !! Code for integrating the blackbody function between two wavenumbers
    !! This is a method that uses a sum convergence.
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

    h = 6.626075540D-34 
    c = 2.99792458D+08
    k = 1.38065812D-23
    c2 = c**2

    do i = 1, nlay1

      if (Tint(i) < 1e-6_dp) then
        bint(i) = 0.0_dp
        cycle
      end if

      do w = 1, 2

        c1 = (h * c) / k
        x = c1 * 100.0_dp * wn_edges(w)/ Tint(i)
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera < 512) then
          itera = 512
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,kind=dp)
          summ = summ + exp(-min(real(j,kind=dp)*x,300.0_dp))*&
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        n2 = 2.0_dp * h * c2
        iB(w) = n2 * (Tint(i)/c1)**4 * summ
      end do
      bint(i) = max(iB(2) - iB(1),0.0_dp)
    end do

  end subroutine BB_integrate

end module ts_short_char_mod_Bezier_bg
