!!!
! Elspeth KH Lee - Aug 2023 : Initial version
!
! lw: Variational Iteration Method four stream variant (alpha-4VIM) - Follows Zhang et al. (2017)
!     Uses (E)AA as an intial guess, then applies VIM to calculate the scattering component
! Pros: Very fast method with LW scattering approximation, no matrix inversions - better than AA alone
! Cons: Still an approximation (though quite good for a non-matrix method)
!!!

module lw_VIM_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these.
  
  !! Optimised quadrature for 1 node - 2 stream (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes - 4 stream (Hogan 2024)
  !integer, parameter :: nmu = 2
  !real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473981_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  !! Gauss–Jacobi-5 quadrature for 3 nodes - 6 stream (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1024922169_dp, 0.4417960320_dp, 0.8633751621_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0437820218_dp, 0.3875796738_dp, 0.5686383044_dp/)

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)

  private :: lw_VIM_fs
  public :: lw_VIM

contains


  subroutine lw_VIM(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nlev) :: be
    real(dp) :: be_int

    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_VIM_fs(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_VIM

  subroutine lw_VIM_fs(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, i, j
    real(dp), dimension(nlay) :: w0, hg, fc
    real(dp), dimension(nlay) :: dtau, beta, epsg, eps, dtau_eg, dtau_e
    real(dp), dimension(nmu, nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nmu, nlay) :: T_eg, T_e, T, cp, cm, wconst, Sp, Sm
    real(dp), dimension(nmu, nmu, nlay) :: Spij, Smij
    real(dp), dimension(nmu, nmu, nlay) :: phip, phim

    real(dp) :: bp, bm, dpp, dm, zepp, zemm, zepm, zemp
    real(dp) :: first, second

    real(dp), dimension(nlay) :: sigma_sq, pmom2, c
    integer, parameter :: nstr = nmu*2
    
    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling

    where (g_in(:) >= 1e-6_dp)
      fc(:) = g_in(:)**(nstr)
      pmom2(:) = g_in(:)**(nstr+1)     
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

    !! Log B with tau function
    !where (dtau(:) < 1.0e-9_dp)
    !  beta(:) = 0.0_dp
    !elsewhere
    beta(:) = -log(be(2:nlev)/be(1:nlay))/dtau(:)
    !end where

    !! modified co-albedo epsilon
    epsg(:) = sqrt((1.0_dp - w0(:))*(1.0_dp - hg(:)*w0(:)))
    !epsg(:) = (1.0_dp - w0(:))

    !! Absorption/modified optical depth for transmission function
    dtau_eg(:) = epsg(:)*dtau(:)

    !! Efficency variables and loop
    do i = 1, nmu
      T_eg(i,:) = exp(-dtau_eg(:)/uarr(i)) ! eg Transmission function
    end do

    !! Start loops to integrate in mu space
    do i = 1, nmu

      !! Begin two-stream loops
      !! Perform downward loop first - also calculate efficency variables
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(i,1) = 0.0_dp
      do k = 1, nlay
        !! Downward AA sweep
        lw_down_g(i,k+1) = lw_down_g(i,k)*T_eg(i,k) + &
          & epsg(k)/(uarr(i)*beta(k) - epsg(k)) * (be(k)*T_eg(i,k) - be(k+1))
      end do

      !! Perform upward loop
      !if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        !lw_up_g(i,nlev) = lw_down_g(i,nlev)*lw_a_surf + (1.0_dp - lw_a_surf)*be_int
      !else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(i,nlev) = lw_down_g(i,nlev) + be_int
      !end if

      do k = nlay, 1, -1
        !! Upward AA sweep
        lw_up_g(i,k) = lw_up_g(i,k+1)*T_eg(i,k) + &
          & epsg(k)/(uarr(i)*beta(k) + epsg(k)) * (be(k) - be(k+1)*T_eg(i,k))
      end do

    end do

    !! If no scattering component in profile, then just find flux and return
    !! no need to perform any scattering calculations
    if (all(w0(:) <= 1.0e-6_dp)) then

      ! Zero the total flux arrays
      flx_up(:) = 0.0_dp
      flx_down(:) = 0.0_dp

      do i = 1, nmu
        !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
        flx_down(:) = flx_down(:) + lw_down_g(i,:) * w(i)
        flx_up(:) = flx_up(:) + lw_up_g(i,:) * w(i)
      end do

      !! The flux is the integrated intensity * pi
      flx_down(:) = pi * flx_down(:)
      flx_up(:) = pi * flx_up(:)

      return

    end if

    !! There is a scattering component, calculate the Sp and Sm

    !! Find Sp and Sm - it's now best to put mu into the inner loop
    ! Sp and Sm defined at lower level edges, zero upper boundary condition
    do k = 1, nlay

      !! Set Sp and Sm to 0
      Sp(:,k) = 0.0_dp
      Sm(:,k) = 0.0_dp

      !! Efficency transmission function calcualtion
      do i = 1, nmu
        T(i,k) = exp(-dtau(k)/uarr(i)) ! regular Transmission function
      end do

      !! co-albedo
      eps(k) = 1.0_dp - w0(k)

      ! Cycle if no scattering component
      if (w0(k) <= 1.0e-6_dp) then
        cycle
      end if

      !! Efficency variables

      !! co-albedo optical depth    
      dtau_e(k) = eps(k)*dtau(k)

      do i = 1, nmu
        T_e(i,k) = exp(-dtau_e(k)/uarr(i)) ! e Transmission function
        cp(i,k) = eps(k)/(uarr(i)*beta(k) + eps(k)) !c+
        cm(i,k) = eps(k)/(-(uarr(i))*beta(k) + eps(k)) !c-
        wconst(i,k) = w0(k)/(real(nmu*2,dp)*(uarr(i))) !constant factor for scattering component
        do j = 1, nmu
          phip(i,j,k) = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(j)   ! phi (net positive mu)
          phim(i,j,k) = 1.0_dp + 3.0_dp*hg(k)*-(uarr(i))*uarr(j) ! phi (net negative mu)
        end do
      end do

      !! Main scattering source function loop
      do i = 1, nmu
        do j = 1, nmu

          !! Note, possible negative sign mistake in Zhang et al. (2017) - zepm and zepp must be positive quantities
          !! To get back the correct expression for the two-stream version
          zepp = -(uarr(i)*uarr(j))/(uarr(i)*eps(k) - uarr(j))
          zemp = (-(uarr(i))*uarr(j))/(-(uarr(i))*eps(k) - uarr(j))
          zepm = -(uarr(i)*-(uarr(j)))/(uarr(i)*eps(k) + uarr(j))
          zemm = (uarr(i)*uarr(j))/(-(uarr(i))*eps(k) + uarr(j))

          first = phip(i,j,k) * zepm * (lw_down_g(j,k) - be(k)*cm(j,k)) * &
            & (1.0_dp - exp(-dtau(k)/zepm))
          second = phim(i,j,k) * zepp * (lw_up_g(j,k+1) - be(k+1)*cp(j,k)) * &
            & (T_e(j,k) - T(i,k))

          Spij(i,j,k) = first + second

          first = phip(i,j,k) * zemp * (lw_up_g(j,k+1) - be(k+1)*cp(j,k)) * &
            & (1.0_dp - exp(-dtau(k)/zemp))
          second = phim(i,j,k) * zemm * (lw_down_g(j,k) - be(k)*cm(j,k)) * &
            & (T_e(j,k) - T(i,k))

          Smij(i,j,k) = first + second

        end do
      end do

      !! Sum up the j Sp/Sm to get the scattering source function
      do i = 1, nmu
        bp = uarr(i)/(uarr(i)*beta(k) + 1.0_dp)
        bm = -(uarr(i))/(-(uarr(i))*beta(k) + 1.0_dp)
        do j = 1, nmu
          Sp(i,k) = Sp(i,k) + &
            & (Spij(i,j,k) - bp*(cp(j,k)*phip(i,j,k) + cm(j,k)*phim(i,j,k))*(be(k+1)*T(i,k) - be(k)))          
          Sm(i,k) = Sm(i,k) + &
            & (Smij(i,j,k) - bm*(cp(j,k)*phim(i,j,k) + cm(j,k)*phip(i,j,k))*(be(k+1) - be(k)*T(i,k)))
        end do
        Sp(i,k) = wconst(i,k)*Sp(i,k)
        Sm(i,k) = wconst(i,k)*Sm(i,k)
      end do

    end do


    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Do final two sweeps including scattering source function - 
    !! Note, don't use AA here, just regular transmission function
    do i = 1, nmu

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(i,1) = 0.0_dp
      do k = 1, nlay

        dm = eps(k)/(-(uarr(i))*beta(k) + 1.0_dp)

        lw_down_g(i,k+1) = lw_down_g(i,k)*T(i,k) + &
          & dm*(be(k+1) - be(k)*T(i,k)) + Sm(i,k)

      end do

      !! Perform upward loop
      !if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        !lw_up_g(i,nlev) = lw_down_g(i,nlev)*lw_a_surf + (1.0_dp - lw_a_surf)*be_int
      !else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(i,nlev) = lw_down_g(i,nlev) + be_int
      !end if

      do k = nlay, 1, -1

        dpp = eps(k)/(uarr(i)*beta(k) + 1.0_dp)

        lw_up_g(i,k) = lw_up_g(i,k+1)*T(i,k) - &
          & dpp*(be(k+1)*T(i,k) - be(k)) + Sp(i,k)

      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(i,:) * w(i)
      flx_up(:) = flx_up(:) + lw_up_g(i,:) * w(i)

    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_VIM_fs

end module lw_VIM_mod