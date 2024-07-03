!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Feautrier
!     Pros: 
!     Cons: 
!!!

module lw_Feautrier_mod
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
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)


  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  private :: lw_Feautrier_method, dtridgl
  public :: lw_Feautrier

contains


  subroutine lw_Feautrier(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_Feautrier_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_net(:))

    !! Net lw fluxes at each level
     lw_up(:) = 0.0_dp
     lw_down(:) = 0.0_dp

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_Feautrier

  subroutine lw_Feautrier_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_net)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_net

    !! Work variables
    integer :: i, k, m
    integer :: niter = 10
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlev) :: w0, hg, Jmean, fd, gd, Sd, flux_tmp_old, flux_tmp, conv_val

    real(dp), dimension(nlev,nmu) :: A, B, C, R, j, h, h_old, l_loc, S
    real(dp) :: dtau_2, d1, d2, f1, f2, f3

    real(dp) :: inv_dtau_min, inv_dtau_min_half

    inv_dtau_min = 1e10_dp
    inv_dtau_min_half = inv_dtau_min * 0.5_dp


    w0(1:nlay) = w_in(:)
    hg(1:nlay) = g_in(:)

    w0(nlev) = 0.0_dp
    hg(nlev) = 0.0_dp

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! Only need the diagonal elements of the matrix, so arrays are nlev by nmu
    A(:,:) = 0.0_dp
    B(:,:) = 0.0_dp
    C(:,:) = 0.0_dp

    !! Initial guess for source function is S = B
    do m = 1, nmu
      S(:,m) = (1.0_dp - w0(:))*be(:)
    end do

    !! Approximate lambda operator
    l_loc(:,:) = 0.0_dp

    !! Find the matrix elements and initial source function
    !! Special upper boundary conditions
    do m = 1, nmu
      f1 = uarr(m)/dtau(1)
      f1 = min(f1,inv_dtau_min)
      A(1,m) = 0.0_dp
      B(1,m) = 1.0_dp + (2.0_dp * f1) * (1.0_dp + f1)
      C(1,m) = -2.0_dp * f1**2
      l_loc(1,m) = l_loc(1,m) + w(m) / (1.0_dp + 2.0_dp*f1 * (1.0_dp + f1))
    end do

    !! Special lower boundary conditions
    do m = 1, nmu
      f1 = uarr(m)/dtau(nlay)
      f1 = min(f1,inv_dtau_min)
      A(nlev,m) = -2.0_dp * f1**2
      B(nlev,m) = 1.0_dp + (2.0_dp * f1) * (1.0_dp + f1)
      C(nlev,m) = 0.0_dp
      l_loc(nlev,m) = l_loc(nlev,m) + w(m) / (1.0_dp + 2.0_dp*f1**2)
    end do

    do k = 2, nlev-1
      dtau_2 = dtau(k)+(dtau(k-1)) 
      do m = 1, nmu

        f1 = 2.0_dp*uarr(m)/dtau_2
        f1 = min(f1,inv_dtau_min_half)
        f2 = uarr(m)/dtau(k)
        f2 = min(f2,inv_dtau_min)
        f3 = uarr(m)/dtau(k-1)
        f3 = min(f3,inv_dtau_min)
        
        A(k,m) = -f1 * f3
        C(k,m) = -f1 * f2
        B(k,m) = 1.0_dp + f1 * (f2 + f3)
        l_loc(k,m) = l_loc(k,m) +  w(m) / (1.0_dp + f1 * (f2 + f3)) 
      end do 
    end do

    !stop

    !! Now we need to converge the source function through
    do i = 1, niter

      !! We now have the matrix coefficents and thermal source function for each level
      !! Solve the tridiagonal system for each angle - updating R to the new source function

      do m = 1, nmu
        R(:,m) = S(:,m)
        R(nlev,m) = R(nlev,m) + (2.0_dp*uarr(m)/dtau(nlay))*be_int
        call dtridgl(nlev,A(:,m),B(:,m),C(:,m),R(:,m),j(:,m))
      end do

      !! Now we use the variable Eddington factor method to iterate to numerical convergence the source function
      !! Find the mean intensity, 1st Eddington coefficent and 2nd Eddington coefficent
      !! and angle integrated source function - assume isotropic source function for now
      Jmean(:) = 0.0_dp
      fd(:) = 0.0_dp
      gd(:) = 0.0_dp
      do m = 1, nmu
        Jmean(:) = Jmean(:) + w(m) * j(:,m)
        fd(:) = fd(:) + w(m) * uarr(m)**2 * j(:,m)
        gd(:) = gd(:) + w(m) * uarr(m) * j(:,m)
      end do
      fd(:) = fd(:)/Jmean(:)
      gd(:) = gd(:)/Jmean(:)

      !! Perform the Variable Eddington Factor method
      do m = 1, nmu
        h(1,m) =-j(1,m)
        h(nlev,m) = 0.0_dp
        do k = 2, nlev - 1
          d1 = (uarr(m)/dtau(k)) * (j(k+1,m) - j(k,m))
          d2 = (uarr(m)/dtau(k-1)) * (j(k,m) - j(k-1,m))
          h(k,m) = -(d1 + d2) * 0.5_dp
        end do
      end do

      flux_tmp_old(:) = flux_tmp(:)
      flux_tmp(:) = 0.0_dp

      do m = 1, nmu
        flux_tmp(:) = flux_tmp(:) - h(:,m)*w(m)*uarr(m)
      end do

      conv_val(:) = abs(1.0_dp - flux_tmp_old(:) / flux_tmp(:))

      if ((all(conv_val(:) < 1e-3_dp))) then
        exit
      end if

      do m = 1, nmu
        S(:,m) = (1.0_dp - w0(:))*be(:) + w0(k)*(0.0_dp + 0.0_dp - l_loc(:,m)*be(:)) &
          & /(1.0_dp - (w0(:) * l_loc(:,m)))
      end do

    end do

    flx_net(:) = 0.0_dp
    do m = 1, nmu
      flx_net(:) = flx_net(:) +  w(m) * uarr(m) * h(:,m) 
      !flx_net(:) = gd(:) * Jmean(:)
    end do
    flx_net(:) = -4.0_dp * pi *  flx_net(:)

    ! Net flux at lower boundary should be internal temperature
    flx_net(nlev) = pi * be_int


    do k = 1, nlev
      !print*, k, j(k,:), Jmean(k), fd(k), gd(k), Sd(k), flx_net(k)
    end do

   !stop

  end subroutine lw_Feautrier_method

  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    integer, parameter :: nmax = 301
    real(dp) :: x
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

end module lw_Feautrier_mod