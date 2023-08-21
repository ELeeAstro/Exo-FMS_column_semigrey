subroutine sw_SDA(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
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
    real(dp), dimension(nlay) :: w0, dtau, hg, T
    real(dp) :: mu1, mu2, mu3, f0
    real(dp) :: om0, om1, om2, om3
    real(dp) :: a0, a1, a2, a3, b0, b1, b2, b3
    real(dp) :: beta, gam, k1, k2, R1, R2, P1, P2, Q1, Q2
    real(dp) :: eta0, eta1, eta2, eta3, del0, del1, del2, del3, delta
    real(dp) :: z1p, z1m, z2p, z2m

    real(dp), dimension(4) :: H1

    !! Calculate dtau in each layer
    dtau(:) = tau_Ve(2:) - tau_Ve(1:nlay)

    !! Delta eddington scaling
    w0(:) = (1.0_dp - g_in(:)**2)*w_in(:)/(1.0_dp-w_in(:)*g_in(:)**2)
    dtau(:) = (1.0_dp-w_in(:)*g_in(:)**2)*dtau(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w0(:) <= 1.0e-3_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    !! Start SDA calculation

    !! First find the Reflection and Transmission coefficents for each layer
    do k = 1, nlay

      !! Transmission to lower level boundary
      T(k) = exp(-tau_Ve(k+1)/mu_z)

      !! Mu moments
      mu1 = mu_z
      mu2 = mu_z**2
      mu3 = mu_z**3

      f0 = 1.0_dp/mu_z

      ! Omega Legendre polynomial 
      om0 = w0(k)
      om1 = 3.0_dp * hg(k)*w0(k) 
      om2 = 5.0_dp * hg(k)**2*w0(k)
      om3 = 7.0_dp * hg(k)**3*w0(k)

      ! Find the a coefficents
      a0 =  1.0_dp - om0
      a1 =  3.0_dp - om1
      a2 =  5.0_dp - om2
      a3 =  7.0_dp - om3

      ! Find the b coefficents
      b0 = 0.25_dp * om0
      b1 = -0.25_dp * om1 * mu1
      b2 = 0.125_dp * om2 * (3.0_dp * mu2 - 1.0_dp)
      b3 = -0.125_dp * om3 * (5.0_dp * mu3 - 3.0d0 * mu1)

      ! Find beta and gamma
      beta = a0*a1 + (4.0_dp/9.0_dp)*a0*a3 + (1.0_dp/9.0_dp)*a2*a3
      gam = (1.0_dp/9.0_dp)*a0*a1*a2*a3

      ! Find k values
      k1 = sqrt(beta + sqrt((beta**2 - 4.0_dp*gam)))/sqrt(2.0_dp)
      k2 = sqrt(beta - sqrt((beta**2 - 4.0_dp*gam)))/sqrt(2.0_dp)

      ! Find R, P and Q coefficents
      R1 = -3.0_dp/2.0_dp * (a0*a1/k1 - k1)/a3
      R2 = -3.0_dp/2.0_dp * (a0*a1/k2 - k2)/a3
      P1 = -a0/k1
      P2 = -a0/k2
      Q1 = 0.5_dp * (a0*a1/k1**2 - 1.0_dp)
      Q2 = 0.5_dp * (a0*a1/k2**2 - 1.0_dp)

      ! Find the delta values
      delta = 9.0_dp*(f0**4 - beta*f0**2 + gam)
      del0 = (a1*b0 - b1*f0)*(a2*a3 - 9.0_dp*f0**2) + 2.0_dp*f0**2*(a3*b2 - 2.0_dp*a3*b0 - 3.0_dp*b3*f0)
      del1 = (a0*b1 - b0*f0)*(a2*a3 - 9.0_dp*f0**2) - 2.0_dp*a0*f0*(a3*b2 - 3.0_dp*b3*f0)
      del2 = (a3*b2 - 3.0_dp*b3*f0)*(a0*a1 - f0**2) - 2.0_dp*a3*f0*(a0*b1 - b0*f0)
      del3 = (a2*b3 - 3.0_dp*b2*f0)*(a0*a1 - f0**2) + f0**2*(6.0_dp*a0*b1 - 4.0_dp*a0*b3 - 6.0_dp*b0*f0)

      ! Find the eta values
      eta0 = del0/delta
      eta1 = del1/delta
      eta2 = del2/delta
      eta3 = del3/delta

      ! Find the Z values
      z1p = 0.5_dp*eta0 + eta1 + 5.0_dp/8.0_dp*eta2
      z1m = 0.5_dp*eta0 - eta1 + 5.0_dp/8.0_dp*eta2
      z2p = -1.0_dp/8.0_dp*eta0 + 5.0_dp/8.0_dp*eta2 + eta3
      z2m = -1.0_dp/8.0_dp*eta0 + 5.0_dp/8.0_dp*eta2 - eta3

      ! Find H vector
      H1(1) = z1m ;  H1(2) = z2m ;  H1(3) = z1p * T(k) ; H1(3) = z2p * T(k)


    end do

  end subroutine sw_SDA