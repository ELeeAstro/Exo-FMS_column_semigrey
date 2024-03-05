!!!
! Elspeth K.H. Lee - 
! sw: 
! lw: 
! Pros:
! Cons:
!!!

module ts_SH_4_mod
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
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Use Two-Term HG function for sw or lw
  logical, parameter :: TTHG_sw = .False.
  logical, parameter :: TTHG_lw = .False.

  public :: ts_SH_4
  private :: sw_SH_4, lw_SH_4, linear_log_interp, bezier_interp

contains

  subroutine ts_SH_4(surf, Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
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
      call sw_SH_4(nlay, nlev, Finc, tau_Ve(:), mu_z(nlev), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
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

    call lw_SH_4(surf, nlay, nlev, pe, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

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

  end subroutine ts_SH_4

  subroutine lw_SH_4(surf, nlay, nlev, pe, be, be_int, tau_IRe, w_in, g_in, lw_a_surf, lw_up, lw_down)
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
    integer, parameter :: nstr = 4
    real(dp), dimension(nstr, nlay) :: a, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface, b_surface_SH4

    real(dp), dimension(nlay) :: beta, gam
    real(dp), dimension(nlay) :: lam1, lam2, exptrm1, exptrm2
    real(dp), dimension(nlay) ::  Q1, Q2, R1, R2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(nlay) :: z1mn_up, z2mn_up, z1pl_up, z2pl_up, z1mn_down, z2mn_down, z1pl_down, z2pl_down
    real(dp), dimension(nlay) :: f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33
    real(dp), dimension(nlev) :: tau_e

    real(dp), parameter :: mu1 = 0.5_dp

    real(dp), dimension(nlay) :: b0
    real(dp), dimension(nlay) :: b1
    real(dp) :: tau_top, b_top

    real(dp), dimension(11,4*nlay) :: Mb
    real(dp), dimension(16,4*nlay) :: Mb_F
    real(dp), dimension(4*nlay) :: B
    real(dp), dimension(nlay) :: X1, X2, X3, X4

    real(dp), dimension(4*nlev, 4*nlay) :: F
    real(dp), dimension(4*nlev) :: G
    real(dp), dimension(4*nlay) :: F_bot
    real(dp) :: G_bot
    
    integer :: info
    integer, dimension(2*nlay) :: ipiv

    integer :: j
    real(dp), dimension(4,4,nlay) :: AA
    real(dp), dimension(4) :: Pubar1
    real(dp), dimension(nlay) :: alpha1, alpha2, beta1, beta2, expo_alp1, expo_alp2, expo_bet1, expo_bet2
    real(dp), dimension(4,nlay) :: exptrm, Aint
    real(dp), dimension(nlay) :: Nint0, Nint1, Nint2, Nint3, expdtau, expo
    real(dp), dimension(nlay) :: multi_scat, intgrl_per_layer
    real(dp) :: u1

    real(dp), dimension(nlev) :: lw_up_g, lw_down_g


    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, g0


    ! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:nlev) - tau_IRe(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling (g**nstream)
    fc(:) = g_in(:)**(nstr)
    pmom2(:) = g_in(:)**(nstr+1)

    ! where (fc(:) /=  pmom2(:))
    !   sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
    !   & ( log(fc(:)**2/pmom2(:)**2) )
    !   c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
    !   fc(:) = c(:)*fc(:)

    !   w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
    !   dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    ! elsewhere
    w0(:) = w_in(:)
    !   fc(:) = 0.0_dp
    ! end where

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
    b_surface_SH4 = (-pi*be(nlev)/4.0_dp)

    surf_reflect = 0.0_dp

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (g0(:) - fc(:)) / (1.0_dp - fc(:))
    w_multi(3,:) = 5.0_dp * (g0(:)**2 - fc(:)) / (1.0_dp - fc(:))
    w_multi(4,:) = 7.0_dp * (g0(:)**3 - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:)
    end do

    !! Find beta and gamma
    beta(:) = a(1,:)*a(2,:) + (4.0_dp/9.0_dp)*a(1,:)*a(4,:) + (1.0_dp/9.0_dp)*a(3,:)*a(4,:)
    gam(:) = (a(1,:)*a(2,:)*a(3,:)*a(4,:))/9.0_dp

    ! Find k values - lambda in Rooney
    lam1(:) = sqrt((beta(:) + sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)
    lam2(:) = sqrt((beta(:) - sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)

    !! Find e values
    exptrm1(:) = min(lam1(:)*dtau(:),35.0_dp)
    exptrm1(:) = exp(-exptrm1(:))
    exptrm2(:) = min(lam2(:)*dtau(:),35.0_dp)
    exptrm2(:) = exp(-exptrm2(:)) 

    R1(:) = -a(1,:)/lam1(:); R2(:) = -a(1,:)/lam2(:)
    Q1(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam1(:)**2) - 1.0_dp)
    Q2(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam2(:)**2) - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a(4,:)) * (a(1,:)*a(2,:)/lam1(:) - lam1(:))
    S2(:) = -3.0_dp/(2.0_dp*a(4,:)) * (a(1,:)*a(2,:)/lam2(:) - lam2(:))

    p1pl(:) = (1.0_dp/2.0_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)  *2.0_dp*pi
    p2pl(:) = (1.0_dp/2.0_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)  *2.0_dp*pi
    q1pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:)) *2.0_dp*pi
    q2pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:)) *2.0_dp*pi
    p1mn(:) = (1.0_dp/2.0_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)  *2.0_dp*pi
    p2mn(:) = (1.0_dp/2.0_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)  *2.0_dp*pi
    q1mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:)) *2.0_dp*pi
    q2mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:)) *2.0_dp*pi

    f00(:)= p1mn(:)*exptrm1(:); f01(:) = p1pl(:)/exptrm1(:); f02(:) = p2mn(:)*exptrm2(:); f03(:) = p2pl(:)/exptrm2(:)
    f10(:) = q1mn(:)*exptrm1(:); f11(:) = q1pl(:)/exptrm1(:); f12(:) = q2mn(:)*exptrm2(:); f13(:) = q2pl(:)/exptrm2(:)
    f20(:) = p1pl(:)*exptrm1(:); f21(:) = p1mn(:)/exptrm1(:); f22(:) = p2pl(:)*exptrm2(:); f23(:) = p2mn(:)/exptrm2(:)
    f30(:) = q1pl(:)*exptrm1(:); f31(:) = q1mn(:)/exptrm1(:); f32(:) = q2pl(:)*exptrm2(:); f33(:) = q2mn(:)/exptrm2(:)

    z1mn_up(:) = (1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp - B1(:)/a(2,:) + B1(:)*dtau(:)/2.0_dp) *2.0_dp*pi 
    z2mn_up(:) = -0.5_dp * (1.0_dp-w0(:)) / (4.0_dp*a(1,:)) * (B0(:) + B1(:)*dtau(:)) *2.0_dp*pi  
    z1pl_up(:) = (1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp + B1(:)/a(2,:) + B1(:)*dtau(:)/2.0_dp) *2.0_dp*pi 
    z2pl_up(:) = -0.5_dp * (1.0_dp-w0(:)) / (4.0_dp*a(1,:)) * (B0(:) + B1(:)*dtau(:)) *2.0_dp*pi  
    z1mn_down(:) = (1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp - B1(:)/a(2,:)) *2.0_dp*pi           
    z2mn_down(:) = -0.5_dp * (1.0_dp-w0(:)) / (4.0_dp*a(1,:)) * (B0(:)) *2.0_dp*pi          
    z1pl_down(:) = (1.0_dp-w0(:))/a(1,:) * (B0(:)/2.0_dp + B1(:)/a(2,:)) *2.0_dp*pi           
    z2pl_down(:) = -0.5_dp * (1.0_dp-w0(:)) / (4.0_dp*a(1,:)) * (B0(:)) *2.0_dp*pi          


    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    !top boundary conditions
    Mb(6,1) = p1mn(1)
    Mb(6,2) = q1pl(1)
    Mb(5,2) = p1pl(1)
    Mb(5,3) = q2mn(1)
    Mb(4,3) = p2mn(1)
    Mb(4,4) = q2pl(1)
    Mb(3,4) = p2pl(1)
    Mb(7,1) = q1mn(1)

    B(1) = b_top - z1mn_down(1)
    B(2) = -b_top/4.0_dp - z2mn_down(1)

    ! bottom boundary conditions
    Mb(6,4*nlay-1) = f22(nlay) - surf_reflect*f02(nlay)
    Mb(6,4*nlay) = f33(nlay) - surf_reflect*f13(nlay)
    Mb(5,4*nlay) = f23(nlay) - surf_reflect*f03(nlay)
    Mb(7,4*nlay-2) = f21(nlay) - surf_reflect*f01(nlay)
    Mb(7,4*nlay-1) = f32(nlay) - surf_reflect*f12(nlay)
    Mb(8,4*nlay-3) = f20(nlay) - surf_reflect*f00(nlay)
    Mb(8,4*nlay-2) = f31(nlay) - surf_reflect*f11(nlay)
    Mb(9,4*nlay-3) = f30(nlay) - surf_reflect*f10(nlay)

    B(4*nlay-1) = b_surface - z1pl_up(nlay) + surf_reflect*z1mn_up(nlay)
    B(4*nlay) = b_surface_SH4 - z2pl_up(nlay) + surf_reflect*z2mn_up(nlay)

    !fill remaining rows of matrix
    do i = 1, nlay-1
      Mb(6,4*i - 1) = f02(i)
      Mb(6,4*i) = f13(i)
      Mb(6,4*i + 1) = -p1pl(i+1)
      Mb(6,4*i + 2) = -q1mn(i+1)
        
      Mb(5,4*i) = f03(i)
      Mb(5,4*i + 1) = -q1mn(i+1)
      Mb(5,4*i + 2) = -p1mn(i+1)
      Mb(5,4*i + 3) = -q2pl(i+1)
        
      Mb(4,4*i + 1) = -p1mn(i+1)
      Mb(4,4*i + 2) = -q1pl(i+1)
      Mb(4,4*i + 3) = -p2pl(i+1)
      Mb(4,4*i + 4) = -q2mn(i+1)
        
      Mb(3,4*i + 2) = -p1pl(i+1)
      Mb(3,4*i + 3) = -q2mn(i+1)
      Mb(3,4*i + 4) = -p2mn(i+1)
        
      Mb(2,4*i + 3) = -p2mn(i+1)
      Mb(2,4*i + 4) = -q2pl(i+1)
        
      Mb(1,4*i + 4) = -p2pl(i+1)
        
      Mb(7,4*i - 2) = f01(i)
      Mb(7,4*i -1) = f12(i)
      Mb(7,4*i) = f23(i)
      Mb(7,4*i + 1) = -q1pl(i+1)
        
      Mb(8,4*i - 3) = f00(i)
      Mb(8,4*i - 2) = f11(i)
      Mb(8,4*i - 1) = f22(i)
      Mb(8,4*i) = f33(i)
        
      Mb(9,4*i - 3) = f10(i)
      Mb(9,4*i - 2) = f21(i)
      Mb(9,4*i - 1) = f32(i)
        
      Mb(10,4*i - 3) = f20(i)
      Mb(10,4*i - 2) = f31(i)
        
      Mb(11,4*i - 3) = f30(i)
    end do 

    do i = 1, nlay-1
      B(4*i - 1) = z1mn_down(i+1) - z1mn_up(i)
      B(4*i) = z2mn_down(i+1) - z2mn_up(i)
      B(4*i + 1) = z1pl_down(i+1)- z1pl_up(i)
      B(4*i + 2) = z2pl_down(i+1) - z2pl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(6,:) = Mb(1,:)
    Mb_F(7,:) = Mb(2,:)
    Mb_F(8,:) = Mb(3,:)  
    Mb_F(9,:) = Mb(4,:)
    Mb_F(10,:) = Mb(5,:)  
    Mb_F(11,:) = Mb(6,:)
    Mb_F(12,:) = Mb(7,:)
    Mb_F(13,:) = Mb(8,:)
    Mb_F(14,:) = Mb(9,:)   
    Mb_F(15,:) = Mb(10,:)    
    Mb_F(16,:) = Mb(11,:)

    call dgbsv(4*nlay, 5, 5, 1, Mb_F, 16, ipiv, B, 4*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

    !! Split B into downward and upward components
    do i = 1, nlay
      X1(i) = B(4*i - 3)
      X2(i) = B(4*i - 2)
      X3(i) = B(4*i - 1)
      X4(i) = B(4*i)
      !print*, i, X1(i), X2(i), X3(i), X4(i)
    end do

    AA(1,1,:) = 1.0_dp; AA(1,2,:) = 1.0_dp; AA(1,3,:) = 1.0_dp;AA(1,4,:) = 1.0_dp
    AA(2,1,:) = R1(:);  AA(2,2,:) = -R1(:); AA(2,3,:) = R2(:); AA(2,4,:) = -R2(:)
    AA(3,1,:) = Q1(:);  AA(3,2,:) =  Q1(:); AA(3,3,:) = Q2(:); AA(3,4,:) =  Q2(:)
    AA(4,1,:) = S1(:);  AA(4,2,:) = -S1(:); AA(4,3,:) = S2(:); AA(4,4,:) = -S2(:)

    !! Now do the source function technique to get up and down fluxes
    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    do m = 1, nmu

      Pubar1(1) = 1.0_dp
      Pubar1(2) = -uarr(m)
      Pubar1(3) = (3.0_dp * -(uarr(m))**2 - 1.0_dp)/2.0_dp
      Pubar1(4) = (5.0_dp * -(uarr(m))**3 - 3.0_dp * -(uarr(m)))/2.0_dp

      u1 = uarr(m)
      alpha1(:) = 1.0_dp/u1 + lam1(:)
      alpha2(:) = 1.0_dp/u1 + lam2(:)
      beta1(:) =  1.0_dp/u1 - lam1(:)
      beta2(:) =  1.0_dp/u1 - lam2(:)
      expo_alp1(:) = min(alpha1(:) * dtau(:),35.0_dp)
      expo_alp2(:) = min(alpha2(:) * dtau(:),35.0_dp)
      expo_bet1(:) = min(beta1(:) * dtau(:) ,35.0_dp)
      expo_bet2(:) = min(beta2(:) * dtau(:) ,35.0_dp)
      exptrm(1,:) = (1.0_dp - exp(-expo_alp1(:))) / alpha1(:) * X1(:)
      exptrm(2,:) = (1.0_dp - exp(-expo_bet1(:))) / beta1(:)  * X2(:)
      exptrm(3,:) = (1.0_dp - exp(-expo_alp2(:))) / alpha2(:) * X3(:)
      exptrm(4,:) = (1.0_dp - exp(-expo_bet2(:))) / beta2(:)  * X4(:)

      Aint(:,:) = 0.0_dp
      do k = 1, nlay
        do j = 1, 4
          do i = 1, 4
            Aint(i, k) = Aint(i,k) + w_multi(j, k) * Pubar1(j) * AA(j, i, k)
          end do
        end do
      end do
      do k = 1, nlay
        do j = 1, 4
          Aint(j,k) = Aint(j,k) * exptrm(j,k)
        end do
      end do


      expdtau(:) = exp(-min(dtau(:)/u1,35.0_dp))
      Nint0(:) = w_multi(1,:) * ((1.0_dp-w0(:)) * u1 / a(1,:) * ( b0(:)*(1.0_dp-expdtau(:)) &
        & + b1(:)*(u1 - (dtau(:)+u1)*expdtau(:))))
      Nint1(:) = w_multi(2,:)*Pubar1(2)* ((1.0_dp-w0(:)) * u1 / a(1,:) * ( b1(:)*(1.0_dp-expdtau(:)) / a(2,:)))
      Nint2(:) = 0.0_dp
      Nint3(:) = 0.0_dp

      multi_scat(:) = Aint(1,:) + Aint(2,:) + Aint(3,:) + Aint(4,:) + Nint0(:) + Nint1(:) + Nint2(:) + Nint3(:)

      expo(:) = min(dtau(:) / uarr(m), 35.0_dp)  
      expdtau(:) = exp(-expo(:))

      intgrl_per_layer(:) = (w0(:) *  multi_scat(:) *2.0_dp*pi &
        & + 2.0_dp*pi*(1.0_dp-w0(:)) * uarr(m) * &
        & (b0(:) * (1.0_dp - expdtau(:)) &
        & + b1(:) * (uarr(m) - (dtau(:) + uarr(m)) * expdtau(:))))

      lw_down_g(1) = 2.0_dp*pi*(1.0_dp - exp(-tau_top/uarr(m)))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = (lw_down_g(k) * exp(-dtau(k)/uarr(m)) + intgrl_per_layer(k) / uarr(m)) 
      end do

      lw_up_g(nlev) = be(nlev) + b1(nlay) * uarr(m)*2.0_dp*pi
      do k = nlay, 1, -1
        lw_up_g(k) = (lw_up_g(k+1) * exp(-dtau(k)/uarr(m)) + intgrl_per_layer(k) / uarr(m)) 
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

    ! print*, lw_up(:)

    ! print*, '-'
    ! print*, lw_down(:)

    !stop

  end subroutine lw_SH_4


  subroutine sw_SH_4(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
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
    integer, parameter :: nstr = 4
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nstr, nlay) :: a, bb, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface, b_top, b_surface_SH4

    real(dp), dimension(nlay) :: beta, gam

    real(dp), dimension(nstr,nlay) :: eta, del
    real(dp), dimension(nlay) :: lam1, lam2, exptrm1, exptrm2
    real(dp), dimension(nlay) :: delta
    real(dp), dimension(nlay) ::  Q1, Q2, R1, R2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(nlay) :: z1mn, z2mn, z1pl, z2pl
    real(dp), dimension(nlay) :: z1mn_up, z2mn_up, z1pl_up, z2pl_up, z1mn_down, z2mn_down, z1pl_down, z2pl_down
    real(dp), dimension(nlay) :: f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33
    real(dp), dimension(nlev) :: expon, tau_e

    real(dp), dimension(11,4*nlay) :: Mb
    real(dp), dimension(16,4*nlay) :: Mb_F
    real(dp), dimension(4*nlay) :: B

    real(dp), dimension(4*nlev, 4*nlay) :: F
    real(dp), dimension(4*nlev) :: G
    real(dp), dimension(4*nlay) :: F_bot
    real(dp) :: G_bot

    real(dp), dimension(4*nlev) :: flux_temp
    real(dp) :: flux_bot, f0


    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, g0
    integer :: info
    integer, dimension(4*nlay) :: ipiv

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

    f0 = 1.0_dp/mu_z

    Pu0(1) = 1.0_dp
    Pu0(2) = -mu_z
    Pu0(3) = (3.0_dp * -(mu_z)**2 - 1.0_dp)/2.0_dp
    Pu0(4) = (5.0_dp * -(mu_z)**3 - 3.0_dp * -(mu_z))/2.0_dp

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (g0(:) - fc(:)) / (1.0_dp - fc(:))
    w_multi(3,:) = 5.0_dp * (g0(:)**2 - fc(:)) / (1.0_dp - fc(:))
    w_multi(4,:) = 7.0_dp * (g0(:)**3 - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:)
      bb(l,:) = ((w0(:) * w_multi(l,:)) * Finc * Pu0(l)) / (4.0_dp*pi)
    end do

    surf_reflect = 0.0_dp

    b_surface = 0.0_dp + surf_reflect*mu_z*Finc*exp(-tau_e(nlev)/mu_z)
    b_surface_SH4 = -(0.0_dp + surf_reflect*mu_z*Finc*exp(-tau_e(nlev)/mu_z))/4.0_dp

    b_top = 0.0_dp

    !! Find beta and gamma
    beta(:) = a(1,:)*a(2,:) + (4.0_dp/9.0_dp)*a(1,:)*a(4,:) + (1.0_dp/9.0_dp)*a(3,:)*a(4,:)
    gam(:) = (a(1,:)*a(2,:)*a(3,:)*a(4,:))/9.0_dp

    ! Find k values - lambda in Rooney
    lam1(:) = sqrt((beta(:) + sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)
    lam2(:) = sqrt((beta(:) - sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)

    !! Find the delta values
    delta(:) = 9.0_dp*(f0**4 - beta(:)*f0**2 + gam(:))
    del(1,:) = (a(2,:)*bb(1,:) - bb(2,:)*f0)*(a(3,:)*a(4,:) - 9.0_dp*f0**2) &
     & + 2.0_dp*f0**2*(a(4,:)*bb(3,:) - 2.0_dp*a(4,:)*bb(1,:) - 3.0_dp*bb(4,:)*f0)
    del(2,:) = (a(1,:)*bb(2,:) - bb(1,:)*f0)*(a(3,:)*a(4,:) - 9.0_dp*f0**2) & 
      & - 2.0_dp*a(1,:)*f0*(a(4,:)*bb(3,:) - 3.0_dp*bb(4,:)*f0)
    del(3,:) = (a(4,:)*bb(3,:) - 3.0_dp*bb(4,:)*f0)*(a(1,:)*a(2,:) & 
      & - f0**2) - 2.0_dp*a(4,:)*f0*(a(1,:)*bb(2,:) - bb(1,:)*f0)
    del(4,:) = (a(3,:)*bb(4,:) - 3.0_dp*bb(3,:)*f0)*(a(1,:)*a(2,:) - f0**2) &
      &  + 2.0_dp*f0**2*(3.0_dp*a(1,:)*bb(2,:) - 2.0_dp*a(1,:)*bb(4,:) - 3.0_dp*bb(1,:)*f0)

    eta(1,:) = del(1,:)/delta(:)
    eta(2,:) = del(2,:)/delta(:)
    eta(3,:) = del(3,:)/delta(:)
    eta(4,:) = del(4,:)/delta(:)


    z1pl(:) = (eta(1,:)/2.0_dp + eta(2,:) + 5*eta(3,:)/8.0_dp)*2.0_dp*pi
    z1mn(:) = (eta(1,:)/2.0_dp - eta(2,:) + 5*eta(3,:)/8.0_dp) *2.0_dp*pi
    z2pl(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8 + eta(4,:))*2.0_dp*pi 
    z2mn(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8 - eta(4,:))*2.0_dp*pi

    !! Find e values
    exptrm1(:) = min(lam1(:)*dtau(:),35.0_dp)
    exptrm1(:) = exp(-exptrm1(:))
    exptrm2(:) = min(lam2(:)*dtau(:),35.0_dp)
    exptrm2(:) = exp(-exptrm2(:))  

    R1(:) = -a(1,:)/lam1(:); R2(:) = -a(1,:)/lam2(:)
    Q1(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam1(:)**2) - 1.0_dp)
    Q2(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam2(:)**2) - 1.0_dp)
    S1(:) = -3.0_dp/(2*a(4,:)) * (a(1,:)*a(2,:)/lam1(:) - lam1(:));
    S2(:) = -3.0_dp/(2.0_dp*a(4,:)) * (a(1,:)*a(2,:)/lam2(:) - lam2(:))

    p1pl(:) = (1.0_dp/2.0_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)  *2.0_dp*pi
    p2pl(:) = (1.0_dp/2.0_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)  *2.0_dp*pi
    q1pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:)) *2.0_dp*pi
    q2pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:)) *2.0_dp*pi
    p1mn(:) = (1.0_dp/2.0_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)  *2.0_dp*pi
    p2mn(:) = (1.0_dp/2.0_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)  *2.0_dp*pi
    q1mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:)) *2.0_dp*pi
    q2mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:)) *2.0_dp*pi

    f00(:)= p1mn(:)*exptrm1(:); f01(:) = p1pl(:)/exptrm1(:); f02(:) = p2mn(:)*exptrm2(:); f03(:) = p2pl(:)/exptrm2(:)
    f10(:) = q1mn(:)*exptrm1(:); f11(:) = q1pl(:)/exptrm1(:); f12(:) = q2mn(:)*exptrm2(:); f13(:) = q2pl(:)/exptrm2(:)
    f20(:) = p1pl(:)*exptrm1(:); f21(:) = p1mn(:)/exptrm1(:); f22(:) = p2pl(:)*exptrm2(:); f23(:) = p2mn(:)/exptrm2(:)
    f30(:) = q1pl(:)*exptrm1(:); f31(:) = q1mn(:)/exptrm1(:); f32(:) = q2pl(:)*exptrm2(:); f33(:) = q2mn(:)/exptrm2(:)

    expon(:) = exp(-tau_e(:)/mu_z)
    z1mn_up(:) = z1mn(:) * expon(2:nlev)
    z2mn_up(:) = z2mn(:) * expon(2:nlev)
    z1pl_up(:) = z1pl(:) * expon(2:nlev)
    z2pl_up(:) = z2pl(:) * expon(2:nlev)
    z1mn_down(:) = z1mn(:) * expon(1:nlay)
    z2mn_down(:) = z2mn(:) * expon(1:nlay)
    z1pl_down(:) = z1pl(:) * expon(1:nlay)
    z2pl_down(:) = z2pl(:) * expon(1:nlay)

    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp
    F_bot(:) = 0.0_dp
    G_bot = 0.0_dp
    F(:,:) = 0.0_dp
    G(:) = 0.0_dp

    !top boundary conditions
    Mb(6,1) = p1mn(1)
    Mb(6,2) = q1pl(1)
    Mb(5,2) = p1pl(1)
    Mb(5,3) = q2mn(1)
    Mb(4,3) = p2mn(1)
    Mb(4,4) = q2pl(1)
    Mb(3,4) = p2pl(1)
    Mb(7,1) = q1mn(1)

    B(1) = b_top - z1mn_down(1)
    B(2) = -b_top/4 - z2mn_down(1)

    ! bottom boundary conditions
    Mb(6,4*nlay-1) = f22(nlay) - surf_reflect*f02(nlay)
    Mb(6,4*nlay) = f33(nlay) - surf_reflect*f13(nlay)
    Mb(5,4*nlay) = f23(nlay) - surf_reflect*f03(nlay)
    Mb(7,4*nlay-2) = f21(nlay) - surf_reflect*f01(nlay)
    Mb(7,4*nlay-1) = f32(nlay) - surf_reflect*f12(nlay)
    Mb(8,4*nlay-3) = f20(nlay) - surf_reflect*f00(nlay)
    Mb(8,4*nlay-2) = f31(nlay) - surf_reflect*f11(nlay)
    Mb(9,4*nlay-3) = f30(nlay) - surf_reflect*f10(nlay)

    B(4*nlay-1) = b_surface - z1pl_up(nlay) + surf_reflect*z1mn_up(nlay)
    B(4*nlay) = b_surface_SH4 - z2pl_up(nlay) + surf_reflect*z2mn_up(nlay)

    !fill remaining rows of matrix
    do i = 1, nlay-1
      Mb(6,4*i - 1) = f02(i)
      Mb(6,4*i) = f13(i)
      Mb(6,4*i + 1) = -p1pl(i+1)
      Mb(6,4*i + 2) = -q1mn(i+1)
        
      Mb(5,4*i) = f03(i)
      Mb(5,4*i + 1) = -q1mn(i+1)
      Mb(5,4*i + 2) = -p1mn(i+1)
      Mb(5,4*i + 3) = -q2pl(i+1)
        
      Mb(4,4*i + 1) = -p1mn(i+1)
      Mb(4,4*i + 2) = -q1pl(i+1)
      Mb(4,4*i + 3) = -p2pl(i+1)
      Mb(4,4*i + 4) = -q2mn(i+1)
        
      Mb(3,4*i + 2) = -p1pl(i+1)
      Mb(3,4*i + 3) = -q2mn(i+1)
      Mb(3,4*i + 4) = -p2mn(i+1)
        
      Mb(2,4*i + 3) = -p2mn(i+1)
      Mb(2,4*i + 4) = -q2pl(i+1)
        
      Mb(1,4*i + 4) = -p2pl(i+1)
        
      Mb(7,4*i - 2) = f01(i)
      Mb(7,4*i -1) = f12(i)
      Mb(7,4*i) = f23(i)
      Mb(7,4*i + 1) = -q1pl(i+1)
        
      Mb(8,4*i - 3) = f00(i)
      Mb(8,4*i - 2) = f11(i)
      Mb(8,4*i - 1) = f22(i)
      Mb(8,4*i) = f33(i)
        
      Mb(9,4*i - 3) = f10(i)
      Mb(9,4*i - 2) = f21(i)
      Mb(9,4*i - 1) = f32(i)
        
      Mb(10,4*i - 3) = f20(i)
      Mb(10,4*i - 2) = f31(i)
        
      Mb(11,4*i - 3) = f30(i)
    end do 

    do i = 1, nlay-1
      B(4*i - 1) = z1mn_down(i+1) - z1mn_up(i)
      B(4*i) = z2mn_down(i+1) - z2mn_up(i)
      B(4*i + 1) = z1pl_down(i+1)- z1pl_up(i)
      B(4*i + 2) = z2pl_down(i+1) - z2pl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(6,:) = Mb(1,:)
    Mb_F(7,:) = Mb(2,:)
    Mb_F(8,:) = Mb(3,:)  
    Mb_F(9,:) = Mb(4,:)
    Mb_F(10,:) = Mb(5,:)  
    Mb_F(11,:) = Mb(6,:)
    Mb_F(12,:) = Mb(7,:)
    Mb_F(13,:) = Mb(8,:)
    Mb_F(14,:) = Mb(9,:)   
    Mb_F(15,:) = Mb(10,:)    
    Mb_F(16,:) = Mb(11,:)


    call dgbsv(4*nlay, 5, 5, 1, Mb_F, 16, ipiv, B, 4*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

   ! flux at bottom of atmosphere
    F_bot(nlay-3) = f20(nlay)
    F_bot(nlay-2) = f21(nlay)
    F_bot(nlay-1) = f22(nlay)
    F_bot(nlay) = f23(nlay)
    G_bot = z1pl_up(nlay)

    F(1,1) = p1mn(1)
    F(1,2) = p1pl(1)
    F(1,3) = p2mn(1)
    F(1,4) = p2pl(1)
    F(2,1) = q1mn(1)
    F(2,2) = q1pl(1)
    F(2,3) = q2mn(1)
    F(2,4) = q2pl(1)
    F(3,1) = p1pl(1)
    F(3,2) = p1mn(1)
    F(3,3) = p2pl(1)
    F(3,4) = p2mn(1)
    F(4,1) = q1pl(1)
    F(4,2) = q1mn(1)
    F(4,2) = q2pl(1)
    F(4,4) = q2mn(1)

    k = 0
    do i = 1, 4*nlay, 4
      F(i+4,i) = f00(k+1)
      F(i+4,i+1) = f01(k+1)
      F(i+4,i+2) = f02(k+1)
      F(i+4,i+3) = f03(k+1)
      F(i+5,i) = f10(k+1)
      F(i+5,i+1) = f11(k+1)
      F(i+5,i+2) = f12(k+1)
      F(i+5,i+3) = f13(k+1)
      F(i+6,i) = f20(k+1)
      F(i+6,i+1) = f21(k+1)
      F(i+6,i+2) = f22(k+1)
      F(i+6,i+3) = f23(k+1)
      F(i+7,i) = f30(k+1)
      F(i+7,i+1) = f31(k+1)
      F(i+7,i+2) = f32(k+1)
      F(i+7,i+3) = f33(k+1)
      k = k + 1
    end do
                             
    G(1) = z1mn_down(1)
    G(2) = z2mn_down(1)
    G(3) = z1pl_down(1)
    G(4) = z2pl_down(1)
    do i = 1, nlay
      G(4*i + 1) = z1mn_up(i)
      G(4*i + 2) = z2mn_up(i)
      G(4*i + 3) = z1pl_up(i)
      G(4*i + 4) = z2pl_up(i)
    end do

    flux_temp(:) = matmul(F(:,:),B(:)) + G(:)
    flux_bot = sum(F_bot(:)*B(:)) + G_bot

    do i = 1, nlev
      sw_down(i) = max(flux_temp(i*4 - 3),0.0_dp) + mu_z * Finc * expon(i)
      sw_up(i) = max(flux_temp(i*4 - 1), 0.0_dp)
    end do

    ! print*, sw_down(:)
    ! print*, sw_up(:)

    ! stop

    !! Now add these values to the diagonal matrix representation
    !! Rooney/PICASO uses a banded matrix structure (which seems like the best idea)
    !! PICASO helpfully figured out the banded matrix structure coefficents, so we follow their indexing directly
    !! The main banded Matix (Mb) is 5x(2*nlay), 2 streams, so 5 terms needed and 2*nlay for the 2 Flux+ terms
    

  end subroutine sw_SH_4

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

end module ts_SH_4_mod
