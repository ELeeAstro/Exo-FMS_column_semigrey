!!!
! Elspeth KH Lee - Nov 2026 : Initial version (shortwave counterpart of lw_sc_scat_ric)
! !!! NOT WELL TESTED YET, USE WITH ABSOLUTE CAUTION !!!
!
! sw: Linear short characteristics, double-Gauss discrete ordinates, delta-M scaled Henyey-Greenstein
!     phase function, level-averaged single scattering properties for the diffuse field - the same
!     discrete problem and Riccati (invariant-imbedding) N x N recursion as lw_sc_scat_ric - with the
!     thermal source q replaced by the single-scattered stellar beam.
!     Only the azimuthally averaged (m = 0) intensity is needed for fluxes.
!
!     Beam source (directional, per stream i, layer k with its own mu0_k as the zcorr beam path):
!       J(t, -/+mu_i) = w_k F_top,k/(4 pi) P(-/+mu_i, -mu0) exp(-t/mu0),
!       P(-mu_i,-mu0) = sum_l (2l+1) chi_l P_l(mu_i) P_l(mu0),  P(+mu_i,-mu0) = sum_l (-1)^l (2l+1) chi_l P_l(mu_i) P_l(mu0)
!     integrated exactly along each stream across the layer (layer optics, not level-averaged):
!       down (to level k+1): bd = J0 P(-mu_i,-mu0) (D/mu_i) exp(-min(1/mu0,1/mu_i) D) em1x(|1/mu0 - 1/mu_i| D)
!       up   (to level k)  : bu = J0 P(+mu_i,-mu0) (D/mu_i) em1x((1/mu0 + 1/mu_i) D)
!     (no resonance at mu_i = mu0). In the lw_sc_scat_ric relations P q_k + Q q_(k+1) -> bd_k (down)
!     and P q_(k+1) + Q q_k -> bu_k (up); the reflection operator recursion is unchanged.
!
!     Boundaries: no diffuse intensity at the top (d_1 = 0); gas giant, black lower boundary
!     (R_K = 0, s_K = 0). Fluxes = diffuse (2 pi sum w mu I) + direct beam F0 mu0 exp(-tau/mu0).
!     The beam follows the zcorr path as sw_AD: mu0_k = mu_z(k+1) and F_top,k = F0 mu_z(nlev)/mu0_k
!     T_beam(k), so the beam scatters exactly the energy it loses.
!
! Notes:
!   - The linear short characteristic for the diffuse source is second order in tau; at GCM layer
!     thicknesses in strongly scattering clouds (and at low sun) errors of a few percent of the incident
!     flux remain. nsub > 1 splits each layer into nsub optically equal sub-layers internally (fluxes
!     returned on the original levels); errors drop roughly as nsub^-2 (nsub^-1 at sharp cloud edges,
!     from the level-averaged w) at nsub times the cost.
!   - The beam phase vectors are normalised so the beam scatters exactly w of its energy (as sw_AD).
!   - NOT energy conserving: the linear short characteristic discretisation of the diffuse field is not
!     flux-conservative, so a purely scattering (w = 1) column shows a spurious net flux divergence of
!     ~3-5% of the incident flux at nsub = 1 (~1% at nsub = 8). sw_AD and the adding methods conserve
!     energy exactly - prefer those inside a GCM; this scheme is for testing.
!
! Pros: any number of streams, exact HG moments to 2 nmu - 1, no particular solution (no mu0 resonance),
!       one N x N recursion step per layer (cheaper than sw_AD's doubling)
! Cons: second order in tau for the diffuse field (needs nsub > 1 in thick clouds); level-averaged w
!!!

module sw_sc_scat_ric_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi

  !! Number of internal optically equal sub-layers per layer (1 = no sub-layering)
  integer, parameter :: nsub = 1

  !! Double-Gauss quadrature: Gauss-Legendre on each hemisphere [0,1], sum(w) = 1, sum(w*uarr) = 1/2
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! 2 nodes per hemisphere - 4 stream
  ! integer, parameter :: nmu = 2
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)

  !! 3 nodes per hemisphere - 6 stream
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.2777777777777778_dp, 0.4444444444444444_dp, 0.2777777777777778_dp/)

  !! 4 nodes per hemisphere - 8 stream
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = &
    & (/0.0694318442029737_dp, 0.3300094782075719_dp, 0.6699905217924281_dp, 0.9305681557970263_dp/)
  real(dp), dimension(nmu), parameter :: w = &
    & (/0.1739274225687269_dp, 0.3260725774312731_dp, 0.3260725774312731_dp, 0.1739274225687269_dp/)

  !! Total number of streams
  integer, parameter :: nstr = 2*nmu

  private :: sw_sc_scat_ric_solve, em1x
  private :: solve_nopiv, lu_factor_nopiv, lu_solve_nopiv, lu_solve_transpose_nopiv
  private :: solve4_nopiv_5, lu_factor4_nopiv, lu_solve4_5_nopiv, lu_solve4_transpose_4_nopiv
  public :: sw_sc_scat_ric

contains

  subroutine sw_sc_scat_ric(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev), intent(in) :: mu_z
    real(dp), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables for the internal sub-layering
    integer :: k, j, ns
    real(dp), dimension(nlay*nsub+1) :: tau_s, sw_up_s, sw_down_s
    real(dp), dimension(nlay*nsub) :: ssa_s, gg_s, mu_s

    !! Shortwave flux calculation
    if ((mu_z(nlev) > 0.0_dp) .and. (Finc > 0.0_dp)) then

      if (nsub == 1) then
        if (mu_z(nlev) == mu_z(1)) then
          mu_s(:) = mu_z(nlev)
        else
          mu_s(:) = mu_z(2:nlev)
        end if
        call sw_sc_scat_ric_solve(nlay, nlev, Finc, mu_z(nlev), mu_s(:), tau_e(:), ssa(:), gg(:), &
          & sw_down(:), sw_up(:))
      else
        !! Split each layer into nsub optically equal sub-layers with the layer properties
        ns = nlay*nsub
        tau_s(1) = tau_e(1)
        do k = 1, nlay
          do j = 1, nsub
            tau_s((k-1)*nsub+j+1) = tau_e(k) + real(j,dp)/real(nsub,dp)*(tau_e(k+1) - tau_e(k))
            ssa_s((k-1)*nsub+j) = ssa(k)
            gg_s((k-1)*nsub+j) = gg(k)
            if (mu_z(nlev) == mu_z(1)) then
              mu_s((k-1)*nsub+j) = mu_z(nlev)
            else
              mu_s((k-1)*nsub+j) = mu_z(k+1)
            end if
          end do
        end do
        call sw_sc_scat_ric_solve(ns, ns+1, Finc, mu_z(nlev), mu_s(:), tau_s(:), ssa_s(:), gg_s(:), &
          & sw_down_s(:), sw_up_s(:))
        do k = 1, nlev
          sw_down(k) = sw_down_s((k-1)*nsub+1)
          sw_up(k) = sw_up_s((k-1)*nsub+1)
        end do
      end if

    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_sc_scat_ric

  subroutine sw_sc_scat_ric_solve(nlay, nlev, F0_in, mu_bot, mu_l, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, mu_bot
    real(dp), dimension(nlay), intent(in) :: mu_l
    real(dp), dimension(nlev), intent(in) :: tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables and arrays
    integer :: k, l, i, j, m, p
    real(dp) :: del, sw, fc, se, so, J0, a0, D, nrm

    real(dp), dimension(nlay) :: dtau, w0, F_top
    real(dp), dimension(nlev) :: cum_trans, T_beam
    real(dp), dimension(0:nstr-1, nlay) :: chi
    real(dp), dimension(nlev) :: ext, sca, ssa_lev
    real(dp), dimension(0:nstr-1, nlev) :: chi_lev
    logical, dimension(nlev) :: noscat

    real(dp), dimension(0:nstr-1, nmu) :: Pl            ! P_l(u_m), upward directions
    real(dp), dimension(0:nstr-1) :: Pmu0               ! P_l(mu0)
    real(dp), dimension(0:nstr-1, nmu, nmu) :: PPw      ! 1/2 (2l+1) P_l(u_i) P_l(u_j) w_j
    real(dp), dimension(nmu, nmu, nlev) :: Fm, Gm       ! forward / backward scattering blocks

    real(dp), dimension(nmu, nlay) :: edel, psi_u, psi_0
    real(dp), dimension(nmu, nlay) :: bd, bu            ! beam source contributions (down, up)
    real(dp), dimension(nmu) :: pdn, pup

    real(dp), dimension(nmu, nmu, nlev) :: Rr, Tt
    real(dp), dimension(nmu, nlev) :: ss, tv
    real(dp), dimension(nmu, nmu) :: Wm, Wi, Vm, Zm, Aa, Bg, Xm, GR, FR
    real(dp), dimension(nmu, nmu+1) :: Ym
    real(dp), dimension(nmu) :: c1, c2, Gs, Fs, e, pu, p0
    real(dp), dimension(nmu) :: dn, up

    !! Delta-M scaled optical depths, albedo and Henyey-Greenstein Legendre moments, f = g**nstr
    do k = 1, nlay
      if (g_in(k) >= 1.0e-6_dp) then
        fc = g_in(k)**nstr
      else
        fc = 0.0_dp
      end if
      w0(k) = min(max(w_in(k), 0.0_dp), 1.0_dp)
      dtau(k) = (1.0_dp - w0(k)*fc)*(tau_in(k+1) - tau_in(k))
      w0(k) = w0(k)*(1.0_dp - fc)/(1.0_dp - fc*w0(k))
      chi(0,k) = 1.0_dp
      do l = 1, nstr-1
        chi(l,k) = (g_in(k)**l - fc)/(1.0_dp - fc)
      end do
    end do

    !! Direct beam transmission to each level along the (zcorr) beam path, scaled tau
    cum_trans(1) = 0.0_dp
    do k = 1, nlay
      cum_trans(k+1) = cum_trans(k) + dtau(k)/mu_l(k)
    end do
    T_beam(:) = exp(-min(cum_trans(:), 200.0_dp))

    !! Beam flux (normal to the beam) at the top of each layer - energy conserving with zcorr
    F_top(:) = F0_in*mu_bot/mu_l(:)*T_beam(1:nlay)

    !! If no scattering at all, return the direct beam only
    if (all(w0(:) <= 1.0e-12_dp)) then
      flx_down(:) = F0_in*mu_bot*T_beam(:)
      flx_up(:) = 0.0_dp
      return
    end if

    !! Layer -> level single scattering properties (identical to lw_sc_scat_ric)
    ext(:) = 0.0_dp
    sca(:) = 0.0_dp
    chi_lev(:,:) = 0.0_dp
    do k = 1, nlay
      sw = w0(k)*dtau(k)
      ext(k) = ext(k) + dtau(k)
      ext(k+1) = ext(k+1) + dtau(k)
      sca(k) = sca(k) + sw
      sca(k+1) = sca(k+1) + sw
      chi_lev(:,k) = chi_lev(:,k) + sw*chi(:,k)
      chi_lev(:,k+1) = chi_lev(:,k+1) + sw*chi(:,k)
    end do
    do k = 1, nlev
      if (ext(k) > 0.0_dp) then
        ssa_lev(k) = sca(k)/ext(k)
      else
        ssa_lev(k) = 0.0_dp
      end if
      if (sca(k) > 0.0_dp) then
        chi_lev(:,k) = chi_lev(:,k)/sca(k)
      else
        chi_lev(:,k) = 0.0_dp
      end if
      chi_lev(0,k) = 1.0_dp
      noscat(k) = (ssa_lev(k) <= 0.0_dp)
    end do

    !! Legendre polynomials at the upward directions (downward: P_l(-u) = (-1)^l P_l(u))
    Pl(0,:) = 1.0_dp
    if (nstr > 1) Pl(1,:) = uarr(:)
    do l = 1, nstr-2
      Pl(l+1,:) = (real(2*l+1,dp)*uarr(:)*Pl(l,:) - real(l,dp)*Pl(l-1,:))/real(l+1,dp)
    end do
    do j = 1, nmu
      do i = 1, nmu
        do l = 0, nstr-1
          PPw(l,i,j) = 0.5_dp*real(2*l+1,dp)*Pl(l,i)*Pl(l,j)*w(j)
        end do
      end do
    end do

    !! Forward (F) and backward (G) scattering blocks: even/odd Legendre split
    !!   F = sum_even + sum_odd,  G = sum_even - sum_odd
    do k = 1, nlev
      if (noscat(k)) then
        Fm(:,:,k) = 0.0_dp
        Gm(:,:,k) = 0.0_dp
        cycle
      end if
      do j = 1, nmu
        do i = 1, nmu
          se = 0.0_dp
          so = 0.0_dp
          do l = 0, nstr-1, 2
            se = se + chi_lev(l,k)*PPw(l,i,j)
          end do
          do l = 1, nstr-1, 2
            so = so + chi_lev(l,k)*PPw(l,i,j)
          end do
          Fm(i,j,k) = ssa_lev(k)*(se + so)
          Gm(i,j,k) = ssa_lev(k)*(se - so)
        end do
      end do
    end do

    !! Linear short characteristic coefficients per layer and stream (as lw_sc_scat_ric)
    do k = 1, nlay
      do m = 1, nmu
        del = dtau(k)/uarr(m)
        edel(m,k) = exp(-del)
        if (del < 1.0e-6_dp) then
          psi_u(m,k) = del*(0.5_dp - del*(1.0_dp/3.0_dp - del*0.125_dp))
          psi_0(m,k) = del*(0.5_dp - del*(1.0_dp/6.0_dp - del/24.0_dp))
        else
          psi_u(m,k) = (1.0_dp - (1.0_dp + del)*edel(m,k))/del
          psi_0(m,k) = (1.0_dp - edel(m,k)) - psi_u(m,k)
        end if
      end do
    end do

    !! Beam source contributions, integrated exactly along each stream with the layer optics
    do k = 1, nlay

      bd(:,k) = 0.0_dp
      bu(:,k) = 0.0_dp
      if (w0(k) <= 0.0_dp .or. dtau(k) <= 0.0_dp) cycle

      !! Beam phase vectors: P(-mu_i,-mu0) (into the downward streams) and P(+mu_i,-mu0) (upward)
      Pmu0(0) = 1.0_dp
      if (nstr > 1) Pmu0(1) = mu_l(k)
      do l = 1, nstr-2
        Pmu0(l+1) = (real(2*l+1,dp)*mu_l(k)*Pmu0(l) - real(l,dp)*Pmu0(l-1))/real(l+1,dp)
      end do
      pdn(:) = 0.0_dp
      pup(:) = 0.0_dp
      do l = 0, nstr-1
        pdn(:) = pdn(:) + real(2*l+1,dp)*chi(l,k)*Pl(l,:)*Pmu0(l)
        pup(:) = pup(:) + (-1.0_dp)**l*real(2*l+1,dp)*chi(l,k)*Pl(l,:)*Pmu0(l)
      end do

      !! Normalise so the beam scatters exactly w0 of its energy
      nrm = 0.5_dp*sum(w(:)*(pdn(:) + pup(:)))
      if (nrm > 0.0_dp) then
        pdn(:) = pdn(:)/nrm
        pup(:) = pup(:)/nrm
      end if

      J0 = w0(k)*F_top(k)/fourpi
      a0 = 1.0_dp/mu_l(k)
      D = dtau(k)
      do i = 1, nmu
        bd(i,k) = J0*pdn(i)*(D/uarr(i))*exp(-min(a0, 1.0_dp/uarr(i))*D)*em1x(abs(a0 - 1.0_dp/uarr(i))*D)
        bu(i,k) = J0*pup(i)*(D/uarr(i))*em1x((a0 + 1.0_dp/uarr(i))*D)
      end do

    end do

    !! ---- Upward Riccati sweep: u_k = R_k d_k + s_k ----
    !! Lower boundary condition - gas giant, black: nothing is reflected back up
    Rr(:,:,nlev) = 0.0_dp
    ss(:,nlev) = 0.0_dp

    do k = nlev-1, 1, -1
      e(:) = edel(:,k)
      pu(:) = psi_u(:,k)
      p0(:) = psi_0(:,k)

      if (noscat(k) .and. noscat(k+1)) then
        !! Clear layer: F = G = 0 -> W = I, Z = E R, X = I, R_k = E R_{k+1} E (diagonal scaling)
        do i = 1, nmu
          c1(i) = bd(i,k)
          c2(i) = e(i)*ss(i,k+1) + bu(i,k)
        end do
        do j = 1, nmu
          do i = 1, nmu
            Rr(i,j,k) = e(i)*Rr(i,j,k+1)*e(j)
            Tt(i,j,k) = 0.0_dp
          end do
          Tt(j,j,k) = e(j)
        end do
        do i = 1, nmu
          sw = 0.0_dp
          do p = 1, nmu
            sw = sw + Rr(i,p,k+1)*c1(p)
          end do
          ss(i,k) = e(i)*sw + c2(i)
          tv(i,k) = c1(i)
        end do
        cycle
      end if

      !! GR = G_{k+1} R_{k+1},  FR = F_{k+1} R_{k+1},  Gs = G_{k+1} s_{k+1},  Fs = F_{k+1} s_{k+1}
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          se = 0.0_dp
          do p = 1, nmu
            sw = sw + Gm(i,p,k+1)*Rr(p,j,k+1)
            se = se + Fm(i,p,k+1)*Rr(p,j,k+1)
          end do
          GR(i,j) = sw
          FR(i,j) = se
        end do
      end do
      do i = 1, nmu
        sw = 0.0_dp
        se = 0.0_dp
        do p = 1, nmu
          sw = sw + Gm(i,p,k+1)*ss(p,k+1)
          se = se + Fm(i,p,k+1)*ss(p,k+1)
        end do
        Gs(i) = sw
        Fs(i) = se
      end do

      !! W, V, c1, c2 (beam sources bd, bu in place of the thermal P q + Q q terms)
      do j = 1, nmu
        do i = 1, nmu
          Wm(i,j) = -p0(i)*(Fm(i,j,k+1) + GR(i,j))
          Vm(i,j) = e(i)*Rr(i,j,k+1) + pu(i)*(FR(i,j) + Gm(i,j,k+1))
        end do
        Wm(j,j) = Wm(j,j) + 1.0_dp
      end do
      do i = 1, nmu
        c1(i) = bd(i,k) + p0(i)*Gs(i)
        c2(i) = e(i)*ss(i,k+1) + pu(i)*Fs(i) + bu(i,k)
      end do

      !! Factor W once. Obtain Z = V W^-1 from W^T Z^T = V^T,
      !! retaining the factors for the later T and t solves.
      if (nmu == 4) then
        call lu_factor4_nopiv(Wm)
      else
        call lu_factor_nopiv(Wm)
      end if
      do j = 1, nmu
        do i = 1, nmu
          Wi(i,j) = Vm(j,i)
        end do
      end do
      if (nmu == 4) then
        call lu_solve4_transpose_4_nopiv(Wm, Wi)
      else
        call lu_solve_transpose_nopiv(nmu, Wm, Wi)
      end if

      !! Transpose Z^T into Z; Aa = E + P F_k, Bg = P G_k
      do j = 1, nmu
        do i = 1, nmu
          Zm(i,j) = Wi(j,i)
          Aa(i,j) = pu(i)*Fm(i,j,k)
          Bg(i,j) = pu(i)*Gm(i,j,k)
        end do
        Aa(j,j) = Aa(j,j) + e(j)
      end do

      !! X = I - Q F_k - Z Bg,  [Y | rhs] = [Q G_k + Z Aa | Z c1 + c2]
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          se = 0.0_dp
          do p = 1, nmu
            sw = sw + Zm(i,p)*Bg(p,j)
            se = se + Zm(i,p)*Aa(p,j)
          end do
          Xm(i,j) = -p0(i)*Fm(i,j,k) - sw
          Ym(i,j) = p0(i)*Gm(i,j,k) + se
        end do
        Xm(j,j) = Xm(j,j) + 1.0_dp
      end do
      do i = 1, nmu
        sw = 0.0_dp
        do p = 1, nmu
          sw = sw + Zm(i,p)*c1(p)
        end do
        Ym(i,nmu+1) = sw + c2(i)
      end do

      !! R_k = X^-1 Y,  s_k = X^-1 (Z c1 + c2)
      if (nmu == 4) then
        call solve4_nopiv_5(Xm, Ym)
      else
        call solve_nopiv(nmu+1, Xm, Ym)
      end if
      Rr(:,:,k) = Ym(:,1:nmu)
      ss(:,k) = Ym(:,nmu+1)

      !! T_k = W^-1 (Aa + Bg R_k),  t_k = W^-1 (Bg s_k + c1)
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          do p = 1, nmu
            sw = sw + Bg(i,p)*Rr(p,j,k)
          end do
          GR(i,j) = Aa(i,j) + sw
        end do
      end do
      do i = 1, nmu
        sw = 0.0_dp
        do p = 1, nmu
          sw = sw + Bg(i,p)*ss(p,k)
        end do
        Gs(i) = sw + c1(i)
      end do
      Ym(:,1:nmu) = GR(:,:)
      Ym(:,nmu+1) = Gs(:)
      if (nmu == 4) then
        call lu_solve4_5_nopiv(Wm, Ym)
      else
        call lu_solve_nopiv(nmu+1, Wm, Ym)
      end if
      Tt(:,:,k) = Ym(:,1:nmu)
      tv(:,k) = Ym(:,nmu+1)

    end do

    !! ---- Downward sweep ----
    !! Top boundary condition - no diffuse intensity downward from the top boundary
    dn(:) = 0.0_dp
    do k = 1, nlev
      if (k > 1) then
        ! d_k = T_{k-1} d_{k-1} + t_{k-1}
        do i = 1, nmu
          sw = tv(i,k-1)
          do p = 1, nmu
            sw = sw + Tt(i,p,k-1)*up(p)     ! up() holds d_{k-1} here
          end do
          c1(i) = sw
        end do
        dn(:) = c1(:)
      end if
      ! u_k = R_k d_k + s_k
      do i = 1, nmu
        sw = ss(i,k)
        do p = 1, nmu
          sw = sw + Rr(i,p,k)*dn(p)
        end do
        c2(i) = sw
      end do
      flx_down(k) = 0.0_dp
      flx_up(k) = 0.0_dp
      do m = 1, nmu
        flx_down(k) = flx_down(k) + dn(m)*w(m)*uarr(m)
        flx_up(k) = flx_up(k) + c2(m)*w(m)*uarr(m)
      end do
      !! Diffuse fluxes plus the direct beam
      flx_down(k) = twopi*flx_down(k) + F0_in*mu_bot*T_beam(k)
      flx_up(k) = twopi*flx_up(k)
      up(:) = dn(:)
    end do

  end subroutine sw_sc_scat_ric_solve

  !! (1 - exp(-x))/x, stable as x -> 0
  elemental function em1x(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: r
    if (x < 1.0e-4_dp) then
      r = 1.0_dp - 0.5_dp*x + x*x/6.0_dp
    else
      r = (1.0_dp - exp(-x))/x
    end if
  end function em1x

  !! Fully unrolled kernels for the default four ordinates per hemisphere.
  !! The generic routines below remain as a fallback if nmu is changed.
  subroutine solve4_nopiv_5(A, B)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: A, B
    real(dp) :: fac

    fac = A(2,1)/A(1,1)
    A(2,2) = A(2,2) - fac*A(1,2)
    A(2,3) = A(2,3) - fac*A(1,3)
    A(2,4) = A(2,4) - fac*A(1,4)
    B(2,1:5) = B(2,1:5) - fac*B(1,1:5)
    fac = A(3,1)/A(1,1)
    A(3,2) = A(3,2) - fac*A(1,2)
    A(3,3) = A(3,3) - fac*A(1,3)
    A(3,4) = A(3,4) - fac*A(1,4)
    B(3,1:5) = B(3,1:5) - fac*B(1,1:5)
    fac = A(4,1)/A(1,1)
    A(4,2) = A(4,2) - fac*A(1,2)
    A(4,3) = A(4,3) - fac*A(1,3)
    A(4,4) = A(4,4) - fac*A(1,4)
    B(4,1:5) = B(4,1:5) - fac*B(1,1:5)

    fac = A(3,2)/A(2,2)
    A(3,3) = A(3,3) - fac*A(2,3)
    A(3,4) = A(3,4) - fac*A(2,4)
    B(3,1:5) = B(3,1:5) - fac*B(2,1:5)
    fac = A(4,2)/A(2,2)
    A(4,3) = A(4,3) - fac*A(2,3)
    A(4,4) = A(4,4) - fac*A(2,4)
    B(4,1:5) = B(4,1:5) - fac*B(2,1:5)

    fac = A(4,3)/A(3,3)
    A(4,4) = A(4,4) - fac*A(3,4)
    B(4,1:5) = B(4,1:5) - fac*B(3,1:5)

    B(4,1:5) = B(4,1:5)/A(4,4)
    B(3,1:5) = B(3,1:5) - A(3,4)*B(4,1:5)
    B(3,1:5) = B(3,1:5)/A(3,3)
    B(2,1:5) = B(2,1:5) - A(2,3)*B(3,1:5)
    B(2,1:5) = B(2,1:5) - A(2,4)*B(4,1:5)
    B(2,1:5) = B(2,1:5)/A(2,2)
    B(1,1:5) = B(1,1:5) - A(1,2)*B(2,1:5)
    B(1,1:5) = B(1,1:5) - A(1,3)*B(3,1:5)
    B(1,1:5) = B(1,1:5) - A(1,4)*B(4,1:5)
    B(1,1:5) = B(1,1:5)/A(1,1)

  end subroutine solve4_nopiv_5

  subroutine lu_factor4_nopiv(A)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: A

    A(2,1) = A(2,1)/A(1,1)
    A(2,2:4) = A(2,2:4) - A(2,1)*A(1,2:4)
    A(3,1) = A(3,1)/A(1,1)
    A(3,2:4) = A(3,2:4) - A(3,1)*A(1,2:4)
    A(4,1) = A(4,1)/A(1,1)
    A(4,2:4) = A(4,2:4) - A(4,1)*A(1,2:4)

    A(3,2) = A(3,2)/A(2,2)
    A(3,3:4) = A(3,3:4) - A(3,2)*A(2,3:4)
    A(4,2) = A(4,2)/A(2,2)
    A(4,3:4) = A(4,3:4) - A(4,2)*A(2,3:4)

    A(4,3) = A(4,3)/A(3,3)
    A(4,4) = A(4,4) - A(4,3)*A(3,4)

  end subroutine lu_factor4_nopiv

  subroutine lu_solve4_5_nopiv(A, B)
    implicit none
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), intent(inout) :: B

    B(2,1:5) = B(2,1:5) - A(2,1)*B(1,1:5)
    B(3,1:5) = B(3,1:5) - A(3,1)*B(1,1:5)
    B(3,1:5) = B(3,1:5) - A(3,2)*B(2,1:5)
    B(4,1:5) = B(4,1:5) - A(4,1)*B(1,1:5)
    B(4,1:5) = B(4,1:5) - A(4,2)*B(2,1:5)
    B(4,1:5) = B(4,1:5) - A(4,3)*B(3,1:5)

    B(4,1:5) = B(4,1:5)/A(4,4)
    B(3,1:5) = B(3,1:5) - A(3,4)*B(4,1:5)
    B(3,1:5) = B(3,1:5)/A(3,3)
    B(2,1:5) = B(2,1:5) - A(2,3)*B(3,1:5)
    B(2,1:5) = B(2,1:5) - A(2,4)*B(4,1:5)
    B(2,1:5) = B(2,1:5)/A(2,2)
    B(1,1:5) = B(1,1:5) - A(1,2)*B(2,1:5)
    B(1,1:5) = B(1,1:5) - A(1,3)*B(3,1:5)
    B(1,1:5) = B(1,1:5) - A(1,4)*B(4,1:5)
    B(1,1:5) = B(1,1:5)/A(1,1)

  end subroutine lu_solve4_5_nopiv

  subroutine lu_solve4_transpose_4_nopiv(A, B)
    implicit none
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), intent(inout) :: B

    B(1,1:4) = B(1,1:4)/A(1,1)
    B(2,1:4) = B(2,1:4) - A(1,2)*B(1,1:4)
    B(2,1:4) = B(2,1:4)/A(2,2)
    B(3,1:4) = B(3,1:4) - A(1,3)*B(1,1:4)
    B(3,1:4) = B(3,1:4) - A(2,3)*B(2,1:4)
    B(3,1:4) = B(3,1:4)/A(3,3)
    B(4,1:4) = B(4,1:4) - A(1,4)*B(1,1:4)
    B(4,1:4) = B(4,1:4) - A(2,4)*B(2,1:4)
    B(4,1:4) = B(4,1:4) - A(3,4)*B(3,1:4)
    B(4,1:4) = B(4,1:4)/A(4,4)

    B(3,1:4) = B(3,1:4) - A(4,3)*B(4,1:4)
    B(2,1:4) = B(2,1:4) - A(3,2)*B(3,1:4)
    B(2,1:4) = B(2,1:4) - A(4,2)*B(4,1:4)
    B(1,1:4) = B(1,1:4) - A(2,1)*B(2,1:4)
    B(1,1:4) = B(1,1:4) - A(3,1)*B(3,1:4)
    B(1,1:4) = B(1,1:4) - A(4,1)*B(4,1:4)

  end subroutine lu_solve4_transpose_4_nopiv

  !! Solve A X = B in place for nr right-hand sides (B overwritten by X), Gaussian elimination,
  !! no pivoting (W and X are diagonally dominant M-matrix Schur complements for p >= 0)
  subroutine solve_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(inout) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do j = 1, nmu
      piv = 1.0_dp/A(j,j)
      do i = j+1, nmu
        fac = A(i,j)*piv
        if (fac /= 0.0_dp) then
          do c = j+1, nmu
            A(i,c) = A(i,c) - fac*A(j,c)
          end do
          do c = 1, nr
            B(i,c) = B(i,c) - fac*B(j,c)
          end do
        end if
      end do
    end do
    do j = nmu, 1, -1
      piv = 1.0_dp/A(j,j)
      do c = 1, nr
        fac = B(j,c)
        do i = j+1, nmu
          fac = fac - A(j,i)*B(i,c)
        end do
        B(j,c) = fac*piv
      end do
    end do

  end subroutine solve_nopiv

  !! In-place no-pivot LU factorization, A = L U. The strict lower triangle
  !! stores the unit-diagonal L multipliers and the upper triangle stores U.
  subroutine lu_factor_nopiv(A)
    implicit none
    real(dp), dimension(nmu, nmu), intent(inout) :: A

    integer :: i, j, c
    real(dp) :: fac, piv

    do j = 1, nmu-1
      piv = 1.0_dp/A(j,j)
      do i = j+1, nmu
        fac = A(i,j)*piv
        A(i,j) = fac
        if (fac /= 0.0_dp) then
          do c = j+1, nmu
            A(i,c) = A(i,c) - fac*A(j,c)
          end do
        end if
      end do
    end do

  end subroutine lu_factor_nopiv

  !! Solve L U X = B using a previously computed no-pivot LU factorization.
  subroutine lu_solve_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(in) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do i = 2, nmu
      do c = 1, nr
        fac = B(i,c)
        do j = 1, i-1
          fac = fac - A(i,j)*B(j,c)
        end do
        B(i,c) = fac
      end do
    end do
    do i = nmu, 1, -1
      piv = 1.0_dp/A(i,i)
      do c = 1, nr
        fac = B(i,c)
        do j = i+1, nmu
          fac = fac - A(i,j)*B(j,c)
        end do
        B(i,c) = fac*piv
      end do
    end do

  end subroutine lu_solve_nopiv

  !! Solve U^T L^T X = B using a previously computed no-pivot LU factorization.
  subroutine lu_solve_transpose_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(in) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do i = 1, nmu
      piv = 1.0_dp/A(i,i)
      do c = 1, nr
        fac = B(i,c)
        do j = 1, i-1
          fac = fac - A(j,i)*B(j,c)
        end do
        B(i,c) = fac*piv
      end do
    end do
    do i = nmu-1, 1, -1
      do c = 1, nr
        fac = B(i,c)
        do j = i+1, nmu
          fac = fac - A(j,i)*B(j,c)
        end do
        B(i,c) = fac
      end do
    end do

  end subroutine lu_solve_transpose_nopiv

end module sw_sc_scat_ric_mod