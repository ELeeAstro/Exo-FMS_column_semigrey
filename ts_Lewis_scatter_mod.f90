!!!
! Elspeth KH Lee - May 2021
! Scattering code of Neil Lewis, following Pierrehumbert (2010)
! Pros: Has scattering
! Cons: Inaccurate at high optical depths
!!!

module ts_Lewis_scatter_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  real(dp), parameter :: gam = 1.0_dp ! associated with closure for integrating over all angles
  real(dp), parameter :: gamprime = 1.0_dp

  public :: ts_Lewis_scatter
  private :: lw_grey_updown, sw_grey_updown, linear_log_interp

contains

  subroutine ts_Lewis_scatter(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, net_F, olr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

    !! Output variables
    real(dp), intent(out) :: olr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlay) :: bl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      call sw_grey_updown(nlay, nlev, F0, tau_Ve(:), sw_a(:), sw_g(:), mu_z, sw_up(:), sw_down(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    bl(:) = sb * Tl(:)**4  ! Integrated planck function flux at levels
    be_int = sb * Tint**4 ! Integrated planck function flux for internal temperature
    call lw_grey_updown(nlay, nlev, bl, be_int, tau_IRe(:), lw_a(:), lw_g(:), mu_z, lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output olr
    olr = lw_up(1)

  end subroutine ts_Lewis_scatter

  subroutine lw_grey_updown(nlay, nlev, bl, be_int, tau_IR, w0, g, mu_z, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_IR
    real(dp), dimension(nlay), intent(in) :: bl, g, w0
    real(dp), intent(in) :: be_int, mu_z

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k
    real(dp), dimension(nlay) :: gh, gam1, gam2, gamb, gamp, gamm
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlev) :: direct

    !! 3/2 closure hemispheric closure condition for g
    gh(:) = g(:) * (3.0_dp/2.0_dp)

    !! Have to compute lw optical depth as per definition in Pierrehumbert (2010)
    do k = 1, nlay
      dtau(k) = -(tau_IR(k+1) - tau_IR(k))
    end do

    !! gammas (as in Pierrehumbert, 2010)
    gam1(:) = gam*(1.0_dp-gh(:)*w0(:)) + gamprime*(1.0_dp-w0(:))
    gam2(:) = gam*(1.0_dp-gh(:)*w0(:)) - gamprime*(1.0_dp-w0(:))
    gamb(:) = gam1(:) - gam2(:)
    gamp(:) = 0.5_dp*w0(:) - gam*w0(:)*gh(:)*1.0_dp
    gamm(:) = 0.5_dp*w0(:) + gam*w0(:)*gh(:)*1.0_dp

    direct(:) = 0.0_dp

    !! Call the two_stream_solver with matrix inversion
    call two_stream_solver(nlay, 1.0_dp, 1.0_dp, be_int, bl, direct, dtau/2.0_dp, &
    & gam1, gam2, gamb, gamp, gamm, lw_up, lw_down)


  end subroutine lw_grey_updown

  subroutine sw_grey_updown(nlay, nlev, Finc, tau_V, w0, g, mu_z, sw_up, sw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_V
    real(dp), dimension(nlay), intent(in) :: w0, g

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down

    !! Work variables
    integer :: k
    real(dp) :: b_int
    real(dp), dimension(nlay) :: gh, gam1, gam2, gamb, gamp, gamm
    real(dp), dimension(nlay) :: bl, dtau
    real(dp), dimension(nlev) :: direct

    gh(:) = g(:) * (3.0_dp/2.0_dp)

    !! Have to compute sw optical depth as per definition in Pierrehumbert (2010)
    do k = 1, nlay
      dtau(k) = -(tau_V(k+1) - tau_V(k))
    end do

    !! gammas (as in Pierrehumbert, 2010)
    gam1(:) = gam*(1.0_dp-gh(:)*w0(:)) + gamprime*(1.0_dp-w0(:))
    gam2(:) = gam*(1.0_dp-gh(:)*w0(:)) - gamprime*(1.0_dp-w0(:))
    gamb(:) = gam1(:) - gam2(:)
    gamp(:) = 0.5_dp*w0(:) - gam*w0(:)*gh(:)*mu_z
    gamm(:) = 0.5_dp*w0(:) + gam*w0(:)*gh(:)*mu_z

    !! Original 1st order method
    direct(1) = mu_z * Finc
    do k = 1, nlay
      direct(k+1) = direct(k) * (2.0_dp*mu_z + dtau(k))/(2.0_dp*mu_z - dtau(k))
      if (direct(k+1) < 0.0_dp) then
        direct(k+1) = 0.0_dp
      end if
    end do

    !! Exponential dependent method
    !direct(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

    !!!! Zero blackbody radiation for shortwave !!!!
    b_int = 0.0_dp
    bl(:) = 0.0_dp

    !! Call the two_stream_solver with matrix inversion
    call two_stream_solver(nlay, 0.0_dp, mu_z, b_int, bl, direct, dtau/2.0_dp, &
    & gam1, gam2, gamb, gamp, gamm, sw_up, sw_down)

    sw_down(:) = sw_down(:) + direct(:)

  end subroutine sw_grey_updown

  subroutine two_stream_solver(nlevels, alpha_g, coszen, pi_B_surf, pi_B,  I_direct, del_tau, gammaone, &
    gammatwo, gammab, gammaplus, gammaminus, I_up, I_down)

  integer, intent(in) :: nlevels
  real(dp), intent(in)              :: alpha_g, coszen,  pi_B_surf  ! should be zero if longwave
  real(dp), intent(in), dimension(nlevels) :: pi_B,  del_tau, gammaone, gammatwo, &
           gammab, gammaplus, gammaminus
           !pi_B zero if SW, Lzeroexp zero if longwave
  real(dp), intent(in), dimension(nlevels+1) :: I_direct


  real(dp), intent(out), dimension(nlevels+1) :: I_up, I_down

  integer :: N  ! N = 2*(nlevels + 1)
  integer, parameter :: MU = 3, ML = 3 ! number of upper and lower diagonals
  real(dp), dimension(2*(nlevels+1), 2*(nlevels+1)) :: A, B
  real(dp), dimension(2*(nlevels+1)) :: I_comb, source_term
  real(dp), dimension(nlevels) :: I_direct_full
  integer :: i, j, lcounter1, lcounter2, IND, M, k, i1, i2
  integer, dimension(2*(nlevels+1)) :: IPVT

  N = 2*(nlevels + 1)

  I_direct_full = (I_direct(1:nlevels) + I_direct(2:nlevels+1)) / 2.0_dp ! move I_direct to full model levels, i.e. I_1 = (I_1/2 + I_3/2) / 2 (approximate, not sure how accurate)


  ! MAKE MATRIX OF ABSORPTION / SCATTERING COEFFICIENTS - band diagonal matrix
  A = 0.0_dp
  ! fill main diagonal
  lcounter1 = 1
  lcounter2 = 1
  do i = 1, N
  if (i.eq.2) then
  A(i,i) = 1 ! upper bc
  elseif (i.eq.N-1) then
  A(i,i) = 1
  else
  if (mod(i,2).eq.1) then
  A(i,i) = -1 + gammaone(lcounter1) / 2.0_dp * del_tau(lcounter1)
  lcounter1 = lcounter1 + 1
  else
  A(i,i) = 1 - gammaone(lcounter2) / 2.0_dp * del_tau(lcounter2)
  lcounter2 = lcounter2 + 1
  endif
  endif
  enddo

  ! fill first upper and lower diagonal
  lcounter1 = 1
  lcounter2 = 1
  do i = 1, N-1
  if (mod(i,2).eq.1) then
  if (i.eq.N-1) then
  A(i,i+1) = - alpha_g ! lower bc
  else
  A(i,i+1) = -gammatwo(lcounter1) / 2.0_dp * del_tau(lcounter1)
  lcounter1 = lcounter1 + 1
  endif
  !else
  if (i.ne.1) then ! for i = 1, A(2,1) = 0 for upper bc
  A(i+1, i) = gammatwo(lcounter2) / 2.0_dp * del_tau(lcounter2)
  lcounter2 = lcounter2 + 1
  endif
  endif
  enddo

  ! fill second upper and lower diagonal
  lcounter1 = 1
  lcounter2 = 1
  do i = 1, N-2
  if (mod(i,2).eq.1) then
  A(i, i+2) = 1.0_dp + gammaone(lcounter1) / 2.0_dp * del_tau(lcounter1)
  lcounter1 = lcounter1 + 1
  else
  A(i+2, i) = -1 - gammaone(lcounter2) / 2.0_dp * del_tau(lcounter2)
  lcounter2 = lcounter2 + 1
  endif
  enddo

  !fill third upper and lower diagonal
  lcounter1 = 1
  lcounter2 = 1
  do i = 1, N-3
  if (mod(i,2).eq.1) then
  A(i, i+3) = - gammatwo(lcounter1) / 2.0_dp * del_tau(lcounter1)
  lcounter1 = lcounter1 + 1
  !else
  A(i+3, i) = gammatwo(lcounter2) / 2.0_dp * del_tau(lcounter2)
  lcounter2 = lcounter2 + 1
  endif
  enddo

  ! make banded matrix B from A
  M = ML + MU + 1
  do j = 1, N
  i1 = max(1,j-MU)
  i2 = min(N,j+ML)
  do i = i1, i2
  k = i - j + M
  B(k,j) = A(i,j)
  enddo
  enddo

  ! Make source term
  source_term = 0.0_dp
  source_term(1) = del_tau(1) * (gammab(1) * pi_B(1) + gammaplus(1) * I_direct_full(1) / coszen)
  lcounter1 = 2
  lcounter2 = 1
  do i = 3, N-2
  if (mod(i,2).eq.1) then
  source_term(i) = del_tau(lcounter1) * (gammab(lcounter1) * pi_B(lcounter1) + &
  gammaplus(lcounter1) * I_direct_full(lcounter1) / coszen)
  lcounter1 = lcounter1 + 1
  else
  source_term(i) = del_tau(lcounter2) * (-1.0_dp*gammab(lcounter2) * pi_B(lcounter2) - &
  gammaminus(lcounter2) * I_direct_full(lcounter2) / coszen)
  lcounter2 = lcounter2 + 1
  endif
  enddo
  source_term(N-1) = alpha_g * I_direct(nlevels+1) + (1.0_dp - alpha_g) * pi_B_surf
  source_term(N) = del_tau(nlevels) * &
    & (-1.0_dp*gammab(nlevels) * pi_B(nlevels) - gammaminus(nlevels) * I_direct_full(nlevels) / coszen)


  ! call LU factorization
  CALL NSBFAC(B,N,N,ML,MU,IPVT,IND)
  ! solve AX = B for X
  if (IND .ne. 0) then
  write(*,*) 'MATRIX IS ILL-CONDITIONED'
  endif
  CALL NSBSLV(B,N,N,ML,MU,IPVT,source_term,I_comb)

  lcounter1 = 1
  lcounter2 = 1
  do i = 1, N
  if (mod(i,2).eq.1) then
  I_up(lcounter1) = I_comb(i)
  lcounter1 = lcounter1 + 1
  else
  I_down(lcounter2) = I_comb(i)
  lcounter2 = lcounter2 + 1
  endif
  enddo


  end subroutine two_stream_solver

  SUBROUTINE NSBFAC(B,LDB,N,ML,MU,IPVT,IND)
  !-------------------------------------------------------------------
  !     LU factorization of a band matrix (non symmetric) with partial
  !     pivoting.
  !     INPUTS:
  !     B  : banded matrix. The correspondance between full matrix
  !     A(i,j) and band matrix B(k,l) is given by following sequence:
  !     m=ml+mu+1
  !     do j=1,n
  !       i1=max(1,j-mu)
  !       i2=min(n,j+ml)
  !       do i=i1,i2
  !         k=i-j+m
  !         b(k,j)=a(i,j)
  !       enddo
  !     enddo
  !     LDB : 1st dimension of B in calling program (ldb.ge.2ml+mu+1)
  !     N   : size of B                             -----------------
  !     ML  : number of lower diagonals
  !     MU  : number of upper diagonals
  !     OUTPUTS:
  !     B   : banded matrix storing the LU elements, the ML first lines
  !           of which are used for pivoting.
  !     IPVT: integer vector of size N storing the pivoting indices.
  !     IND : flag = 0,  B is non singular
  !                = k,  B may be singular
  !-------------------------------------------------------------------
  INTEGER :: N, LDB, ML, MU, IND, OUT
  INTEGER :: I, M, J0, J1, JZ, I0, J, JU, K, KP1, L, LM, MM, n_m_one
  REAL(dp) :: B(LDB,N),T
  INTEGER IPVT(*)
  M=ML+MU+1
  IND=0

  J0=MU+2
  J1=MIN(N,M)-1
  IF(J1.GE.J0) THEN
  DO JZ=J0,J1
  I0=M+1-JZ
  DO I=I0,ML
  B(I,JZ)=0.0_dp
  ENDDO
  ENDDO
  ENDIF
  JZ=J1
  JU=0

  n_m_one=N-1
  IF(n_m_one.GE.1) THEN
  DO K=1,n_m_one
  KP1=K+1
  JZ=JZ+1
  IF(JZ.LE.N.AND.ML.GE.1) THEN
  DO I=1,ML
  B(I,JZ)=0.0_dp
  ENDDO
  ENDIF

  LM=MIN(ML,N-K)
  call IAMAX(B(M,K),LM+1, OUT)
  L=OUT+M-1
  IPVT(K)=L+K-M
  IF(B(L,K).EQ.0.0_dp) GO TO 10
  IF(L.NE.M) THEN
  T=B(L,K)
  B(L,K)=B(M,K)
  B(M,K)=T
  ENDIF
  T=-1.0_dp/B(M,K)
  CALL SCALE(LM,T,B(M+1,K))

  JU=MIN(MAX(JU,MU+IPVT(K)),N)
  MM=M
  IF(JU.GE.KP1) THEN
  DO J=KP1,JU
  L=L-1
  MM=MM-1
  T=B(L,J)
  IF(L.NE.MM) THEN
  B(L,J)=B(MM,J)
  B(MM,J)=T
  ENDIF
  CALL DAXPY(LM,T,B(M+1,K),B(MM+1,J))
  ENDDO
  ENDIF
  GO TO 20
  10     ind=k
  20     CONTINUE
  ENDDO
  ENDIF
  IPVT(N)=N
  IF(B(M,N).EQ.0.0_dp) IND=N
  END SUBROUTINE NSBFAC

  SUBROUTINE DAXPY(N,A,X,Y)
  INTEGER N, I
  REAL(dp) X(*),Y(*),A
  DO I=1,N
  Y(I)=Y(I)+A*X(I)
  ENDDO
  END SUBROUTINE DAXPY

  SUBROUTINE IAMAX(A,N,FINAL)
  INTEGER N, I, FINAL
  REAL(dp) A(*),T
  T=0.0_dp
  DO I=1,N
  IF(ABS(A(I)).GT.T) THEN
  t=abs(A(I))
  FINAL=I
  ENDIF
  enddo
  END SUBROUTINE IAMAX

  SUBROUTINE SCALE(N,T,A)
  INTEGER N, I
  REAL(dp) A(*),T
  DO I=1,N
  A(I)=T*A(I)
  ENDDO
  END SUBROUTINE SCALE

  SUBROUTINE NSBSLV(A,LDA,N,ML,MU,IPVT,b,x)
  !--------------------------------------------------------------------
  !     Solve banded linear system Ax = b
  !     INPUTS:
  !     A   : banded matrix as output of LU factorization by NSBFAC
  !           (see storing mode in NSBFAC subroutine).
  !     LDA : 1st dimension of A in calling program (lda.ge.2ml+mu+1)
  !     N   : order of A                             -----------------
  !     ML  : number of lower diagonals
  !     MU  : number of upper diagonals
  !     IPVT: integer vector of size N storing the pivoting indices
  !           as output of NSBFAC.
  !     b   : second member vector
  !     OUTPUT:
  !     x   : solution vector
  !---------------------------------------------------------------------

  INTEGER :: N, ML, MU, LDA
  REAL(dp) :: A(LDA,*),B(*),X(*),T
  INTEGER IPVT(*)
  INTEGER :: I, K, L, M, LM, LA, LB, n_m_one
  DO I=1,N
  X(I)=B(I)
  ENDDO
  M=ML+MU+1
  n_m_one=N-1
  !     solve L*y = b
  IF(ML.NE.0.AND.n_m_one.GE.1) THEN
  DO K=1,n_m_one
  LM=MIN(ML,N-K)
  L=IPVT(K)
  T=X(L)
  IF(L.NE.K) THEN
  X(L)=X(K)
  X(K)=T
  ENDIF
  CALL DAXPY(LM,T,A(M+1,K),X(K+1))
  ENDDO
  ENDIF
  !     solve U*y = x
  DO K=N,1,-1
  X(K)=X(K)/A(M,K)
  LM=MIN(K,M)-1
  LA=M-LM
  LB=K-LM
  T=-X(K)
  CALL DAXPY(LM,T,A(LA,K),X(LB))
  ENDDO
  END SUBROUTINE NSBSLV

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

end module ts_Lewis_scatter_mod
