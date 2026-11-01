!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!
! lw: Four-stream discrete ordinates adding method (4DDA; Zhang et al. 2013, J. Atmos. Sci.;
!     infrared adding following the Zhang group's four-stream infrared work), thermal emission only.
!     Layer:   analytic 4-stream (2 double-Gauss streams per hemisphere) discrete ordinates
!              solution of each homogeneous layer. With S = I+ + I-, D = I+ - I-, the system
!              reduces to d2S/dtau2 = (alpha+beta)(alpha-beta) S, a 2x2 eigenproblem whose
!              eigenvalues k1^2, k2^2 are the roots of a quadratic (Zhang et al.'s biquadratic).
!              The Planck function is linear in tau inside the layer: I+/- = B(tau) +/- d B1,
!              with (E - w (Ws - Wd)) d = mu, so the thermal particular solution is exact.
!              Only decaying exponentials are used (no overflow for any dtau).
!     Adding:  each layer gives 2x2 reflection R and transmission T (equal from both sides,
!              homogeneous layer) and 2-vector thermal sources (Eu emitted up at the top, Ed
!              down at the bottom); composites from the top (Rtop, Sdn) and from the bottom
!              (Rbot, Sup, starting from the gas giant boundary I_up = I_down + B(Tint)) are
!              combined at every level - the same structure as sw_SDA and lw_AD.
!
! Notes:
!   - Henyey-Greenstein phase function with delta-M scaling, f = g**4 (chi_(2 nmu), the first
!     moment the 4-stream expansion cannot hold; the retained moments chi_0..chi_3 are exact).
!   - Double-Gauss quadrature (Gauss-Legendre on [0,1], mu = 1/2 -/+ 1/(2 sqrt(3)), w = 1/2),
!     as Zhang et al.; the streams are the scattering streams, so flux = 2 pi sum(w mu I).
!     With 2 nodes per hemisphere sum_j w_j (P(mu_i,mu_j) + P(mu_i,-mu_j))/2 = 1 exactly,
!     so energy is conserved without renormalising the phase matrix and an isothermal
!     layer gives I = B exactly.
!   - w is capped at w_max (the k2 -> 0 conservative limit is singular in the eigen-solution).
!
! Pros: four-stream accuracy for multiple scattering at two-stream-like cost, no iterations,
!       only 4x4 (layer) and 2x2 (adding) solves
! Cons: four streams only - for very strong forward scattering lw_AD with more streams is better
!!!

module lw_DDA_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of streams per hemisphere (4-stream method) and retained phase function moments
  integer, parameter :: nmu = 2
  integer, parameter :: nmom = 2*nmu

  !! Double-Gauss nodes and weights on [0,1] (sum(w) = 1, sum(w*mu) = 1/2)
  real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)

  !! Single scattering albedo cap and minimum layer optical depth
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-6_dp
  real(dp), parameter :: dtau_min = 1.0e-10_dp

  private :: lw_DDA_method, layer_RT, solve4, inv2
  public :: lw_DDA

contains

  subroutine lw_DDA(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: Te, be
    real(dp) :: be_int

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_DDA_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_DDA

  subroutine lw_DDA_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, l, i, j
    real(dp), dimension(nmom, nmu) :: Pl
    real(dp), dimension(nmom) :: chi, c2l, sgn
    real(dp), dimension(nlay) :: dtau, w0
    real(dp), dimension(nmu,nmu) :: Ws, Wd, Imat, Mi
    real(dp), dimension(nmu,nmu,nlay) :: R, T
    real(dp), dimension(nmu,nlay) :: Eu, Ed
    real(dp), dimension(nmu,nmu,nlev) :: Rtop, Rbot
    real(dp), dimension(nmu,nlev) :: Sdn, Sup
    real(dp), dimension(nmu) :: Iup, Idn
    real(dp) :: f

    !! Legendre polynomials P_l(mu_i), l = 0..3 (stored at l+1)
    Pl(1,:) = 1.0_dp
    Pl(2,:) = uarr(:)
    Pl(3,:) = 0.5_dp*(3.0_dp*uarr(:)**2 - 1.0_dp)
    Pl(4,:) = 0.5_dp*(5.0_dp*uarr(:)**3 - 3.0_dp*uarr(:))

    !! (2l+1) and (-1)^l factors of the phase function expansion
    do l = 0, nmom-1
      c2l(l+1) = real(2*l+1,dp)
      sgn(l+1) = (-1.0_dp)**l
    end do

    Imat(:,:) = 0.0_dp
    do i = 1, nmu
      Imat(i,i) = 1.0_dp
    end do

    !! Layer optical depths and single scattering albedo
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), w_max)

    !! Reflection, transmission and thermal sources of each layer
    do k = 1, nlay

      !! Henyey-Greenstein moments, delta-M scaled with the first truncated moment chi_4 = g^4
      do l = 0, nmom-1
        chi(l+1) = g_in(k)**l
      end do
      f = g_in(k)**nmom
      if (w0(k) > 0.0_dp .and. f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
        chi(:) = (chi(:) - f)/(1.0_dp - f)
        dtau(k) = (1.0_dp - w0(k)*f)*dtau(k)
        w0(k) = w0(k)*(1.0_dp - f)/(1.0_dp - w0(k)*f)
      end if

      !! Quadrature-weighted phase matrices, normalised so (1/2) int P dmu' = 1:
      !! Ws(i,j) = w_j P(mu_i, mu_j)/2 (same hemisphere), Wd(i,j) = w_j P(mu_i,-mu_j)/2 (opposite)
      do j = 1, nmu
        do i = 1, nmu
          Ws(i,j) = 0.5_dp*w(j)*sum(c2l(:)*chi(:)*Pl(:,i)*Pl(:,j))
          Wd(i,j) = 0.5_dp*w(j)*sum(sgn(:)*c2l(:)*chi(:)*Pl(:,i)*Pl(:,j))
        end do
      end do

      call layer_RT(dtau(k), w0(k), be(k), be(k+1), Ws, Wd, Imat, R(:,:,k), T(:,:,k), Eu(:,k), Ed(:,k))

    end do

    !! Composites from the top: Rtop = reflection of the slab above a level (for upward light),
    !! Sdn = downward intensity emitted by the slab above - no incident intensity at the top
    Rtop(:,:,1) = 0.0_dp
    Sdn(:,1) = 0.0_dp
    do k = 1, nlay
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rtop(:,:,k), R(:,:,k)))
      Rtop(:,:,k+1) = R(:,:,k) + matmul(T(:,:,k), matmul(Mi(:,:), matmul(Rtop(:,:,k), T(:,:,k))))
      Sdn(:,k+1) = Ed(:,k) + matmul(T(:,:,k), matmul(Mi(:,:), Sdn(:,k) + matmul(Rtop(:,:,k), Eu(:,k))))
    end do

    !! Composites from the bottom: Rbot = reflection of the slab below a level (for downward light),
    !! Sup = upward intensity emitted by the slab below.
    !! Lower boundary - gas giant: I_up = I_down + B(Tint) in every stream
    Rbot(:,:,nlev) = Imat(:,:)
    Sup(:,nlev) = be_int
    do k = nlay, 1, -1
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rbot(:,:,k+1), R(:,:,k)))
      Rbot(:,:,k) = R(:,:,k) + matmul(T(:,:,k), matmul(Mi(:,:), matmul(Rbot(:,:,k+1), T(:,:,k))))
      Sup(:,k) = Eu(:,k) + matmul(T(:,:,k), matmul(Mi(:,:), Sup(:,k+1) + matmul(Rbot(:,:,k+1), Ed(:,k))))
    end do

    !! Intensities at each level from the two composites:
    !! I_dn = Sdn + Rtop I_up,  I_up = Sup + Rbot I_dn
    do k = 1, nlev
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rtop(:,:,k), Rbot(:,:,k)))
      Idn(:) = matmul(Mi(:,:), Sdn(:,k) + matmul(Rtop(:,:,k), Sup(:,k)))
      Iup(:) = Sup(:,k) + matmul(Rbot(:,:,k), Idn(:))
      flx_down(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Idn(:))
      flx_up(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Iup(:))
    end do

  end subroutine lw_DDA_method

  !! Analytic 4-stream discrete ordinates solution of one homogeneous layer:
  !! R, T (2x2, incident from either side) and thermal sources Eu (up, top), Ed (down, bottom)
  subroutine layer_RT(dtau, w0, B_top, B_bot, Ws, Wd, Imat, R, T, Eu, Ed)
    implicit none

    real(dp), intent(in) :: dtau, w0, B_top, B_bot
    real(dp), dimension(nmu,nmu), intent(in) :: Ws, Wd, Imat
    real(dp), dimension(nmu,nmu), intent(out) :: R, T
    real(dp), dimension(nmu), intent(out) :: Eu, Ed

    integer :: i, m
    real(dp), dimension(nmu,nmu) :: al, bt, amb, apb, Q, Gp, Gm
    real(dp), dimension(nmu) :: kk, ex, X, Y, v1, v2, d, Ipm0, Ipp0, IppD, IpmD
    real(dp), dimension(2*nmu,2*nmu) :: L
    real(dp), dimension(2*nmu,3) :: rhs
    real(dp), dimension(nmu,3) :: Up, Dn
    real(dp) :: tr, det, disc, lam1, lam, B1

    !! Optically negligible layer - transparent
    if (dtau <= dtau_min) then
      R(:,:) = 0.0_dp
      T(:,:) = Imat(:,:)
      Eu(:) = 0.0_dp
      Ed(:) = 0.0_dp
      return
    end if

    !! mu dI+/dtau = alpha I+ - beta I-,  dI-/dtau = beta I+ - alpha I-  (per unit tau, divided by mu)
    do i = 1, nmu
      al(i,:) = (Imat(i,:) - w0*Ws(i,:))/uarr(i)
      bt(i,:) = w0*Wd(i,:)/uarr(i)
    end do
    amb(:,:) = al(:,:) - bt(:,:)
    apb(:,:) = al(:,:) + bt(:,:)
    Q(:,:) = matmul(apb(:,:), amb(:,:))

    !! Eigenvalues k^2 of Q (quadratic in k^2); smaller root via det/lam1 to avoid cancellation
    tr = Q(1,1) + Q(2,2)
    det = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1)
    disc = sqrt(max(tr**2 - 4.0_dp*det, 0.0_dp))
    lam1 = 0.5_dp*(tr + disc)
    kk(1) = sqrt(lam1)
    kk(2) = sqrt(max(det/lam1, tiny(1.0_dp)))

    !! Eigenvectors X of Q, then the decaying (e^{-k tau}) modes:
    !! I+ = Gp e^{-k tau}, I- = Gm e^{-k tau}, with Gp = (X - Y)/2, Gm = (X + Y)/2, Y = (alpha-beta) X/k
    !! The growing mode is the mirror image (I+ = Gm, I- = Gp), written as e^{-k (dtau - tau)}
    do m = 1, nmu
      lam = kk(m)**2
      v1(:) = (/Q(1,2), lam - Q(1,1)/)
      v2(:) = (/lam - Q(2,2), Q(2,1)/)
      if (sum(v1**2) >= sum(v2**2)) then
        X(:) = v1(:)/sqrt(sum(v1**2))
      else
        X(:) = v2(:)/sqrt(sum(v2**2))
      end if
      Y(:) = matmul(amb(:,:), X(:))/kk(m)
      Gp(:,m) = 0.5_dp*(X(:) - Y(:))
      Gm(:,m) = 0.5_dp*(X(:) + Y(:))
      ex(m) = exp(-kk(m)*dtau)
    end do

    !! Thermal particular solution, B(tau) = B_top + B1 tau: I+/- = B(tau) +/- d B1
    B1 = (B_bot - B_top)/dtau
    d(:) = matmul(inv2(Imat(:,:) - w0*(Ws(:,:) - Wd(:,:))), uarr(:))
    Ipp0(:) = B_top + d(:)*B1   ! I+ at the top
    Ipm0(:) = B_top - d(:)*B1   ! I- at the top
    IppD(:) = B_bot + d(:)*B1   ! I+ at the bottom
    IpmD(:) = B_bot - d(:)*B1   ! I- at the bottom

    !! Boundary value problem for the coefficients (A1, A2, C1, C2):
    !! I+ = sum_m [A_m Gp_m e^{-k_m tau} + C_m Gm_m e^{-k_m (dtau - tau)}] + Ip+
    !! I- = sum_m [A_m Gm_m e^{-k_m tau} + C_m Gp_m e^{-k_m (dtau - tau)}] + Ip-
    !! Rows 1-2: I-(0) given, rows 3-4: I+(dtau) given
    do m = 1, nmu
      L(1:nmu,m) = Gm(:,m)
      L(1:nmu,nmu+m) = Gp(:,m)*ex(m)
      L(nmu+1:2*nmu,m) = Gp(:,m)*ex(m)
      L(nmu+1:2*nmu,nmu+m) = Gm(:,m)
    end do

    !! RHS 1,2: unit incident intensity in stream 1,2 at the top (no source) -> R, T columns
    !! RHS 3: no incident intensity, thermal source only -> Eu, Ed
    rhs(:,:) = 0.0_dp
    rhs(1,1) = 1.0_dp
    rhs(2,2) = 1.0_dp
    rhs(1:nmu,3) = -Ipm0(:)
    rhs(nmu+1:2*nmu,3) = -IppD(:)

    call solve4(L, rhs)

    !! Outgoing intensities: I+ at the top and I- at the bottom (homogeneous part)
    do m = 1, 3
      Up(:,m) = matmul(Gp(:,:), rhs(1:nmu,m)) + matmul(Gm(:,:), ex(:)*rhs(nmu+1:2*nmu,m))
      Dn(:,m) = matmul(Gm(:,:), ex(:)*rhs(1:nmu,m)) + matmul(Gp(:,:), rhs(nmu+1:2*nmu,m))
    end do

    R(:,1:nmu) = Up(:,1:nmu)
    T(:,1:nmu) = Dn(:,1:nmu)
    Eu(:) = Up(:,3) + Ipp0(:)
    Ed(:) = Dn(:,3) + IpmD(:)

  end subroutine layer_RT

  !! Solve the 4x4 system L x = b for 3 right hand sides (Gaussian elimination, partial pivoting)
  !! On exit b holds the solutions
  subroutine solve4(L, b)
    implicit none

    real(dp), dimension(2*nmu,2*nmu), intent(inout) :: L
    real(dp), dimension(2*nmu,3), intent(inout) :: b

    integer :: n, i, j, p
    real(dp) :: fac
    real(dp), dimension(2*nmu) :: rowL
    real(dp), dimension(3) :: rowb

    n = 2*nmu

    do j = 1, n
      p = j - 1 + maxloc(abs(L(j:n,j)), dim=1)
      if (p /= j) then
        rowL(:) = L(j,:) ; L(j,:) = L(p,:) ; L(p,:) = rowL(:)
        rowb(:) = b(j,:) ; b(j,:) = b(p,:) ; b(p,:) = rowb(:)
      end if
      do i = j+1, n
        fac = L(i,j)/L(j,j)
        L(i,j:n) = L(i,j:n) - fac*L(j,j:n)
        b(i,:) = b(i,:) - fac*b(j,:)
      end do
    end do

    do i = n, 1, -1
      b(i,:) = (b(i,:) - matmul(L(i,i+1:n), b(i+1:n,:)))/L(i,i)
    end do

  end subroutine solve4

  !! Inverse of a 2x2 matrix
  pure function inv2(A) result(B)
    implicit none

    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2,2) :: B
    real(dp) :: det

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    B(1,1) = A(2,2)/det
    B(1,2) = -A(1,2)/det
    B(2,1) = -A(2,1)/det
    B(2,2) = A(1,1)/det

  end function inv2

end module lw_DDA_mod