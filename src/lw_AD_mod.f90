!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Adding-Doubling method
!!!

module lw_AD_mod
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
  !! For HJ's we actually care about stratospheric heating rates, so Gauss–Laguerre quadrature is generally best for 4+ stream, 
  !! even for cloudy regions (Hogan 2024).
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Laguerre quadrature for 2 nodes (Hogan 2024)
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.1813898346_dp, 0.7461018061_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.1464466094_dp, 0.8535533906_dp/)

  !! Gauss–Laguerre quadrature for 3 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.0430681066_dp, 0.3175435896_dp, 0.8122985952_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0103892565_dp, 0.2785177336_dp, 0.7110930099_dp/)

  !! Gauss–Laguerre quadrature for 4 nodes (Hogan 2024)
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.0091177205_dp, 0.1034869099_dp, 0.4177464746_dp, 0.8510589811_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0005392947_dp, 0.0388879085_dp, 0.3574186924_dp, 0.6031541043_dp /)

  !! Optical depth limit on optical depth of layers
  real(dp), parameter :: dtau_lim = 0.01_dp

  private :: lw_doubling_adding, matinv2, matinv4, ludcmp, lubksb, inv_LU
  public :: lw_AD

contains


  subroutine lw_AD(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, a_surf, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: a_surf
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
    call lw_adding_doubling(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AD

  subroutine lw_adding_doubling(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m
    real(dp), dimension(nlay) :: bm, bl
    real(dp), dimension(nlay) :: dtau, w0, hg

    real(dp), dimension(nmu, nmu) :: G
    real(dp), dimension(nmu, nmu) :: I


    real(dp), dimension(nlev, nmu) :: sp, sm, sp_add, sm_add

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    w0(:) = w_in(:)
    hg(:) = g_in(:)


    do k = 1, nlay
      if (dtau(k) <= 1.0e-6_dp) then
        ! For low optical depths use the isothermal approimation
        bm(k) = 0.5_dp*(be(k+1) + be(k))
        bl(k) = 0.0_dp
      else
        bm(k) = 0.5_dp*(be(k+1) + be(k))
        bl(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
      endif
    end do

    ! TOA Reflection boundary
    R(1,:,:) = 
    ! TOA Transmission boundary
    T(1,:,:) = 
    ! TOA emissivity and slope
    y(1,:) =
    z(1,:) = 

    sp(1,:) =
    sm(1,:) =  

    !! Identity matrix
    I(:,:) = 0.0_dp

    !! Begin loop from top of atmosphere to bottom to find transmision and emission and each layer boundary
    !! Perform doubling scheme as required when dtau > dtau_lim
    do k = 1, nlay

      !if (dtau(k) <= dtau_lim*uarr(1)) then
        !! Safe to perform the adding method directly

        gt = dtau(k)

        !! Calculate Gamma 
        G(:,:) = I(:,:) - matmul(R(k,:,:),R(k,:,:))
        G(:,:) = matinv2(G(:,:)) ! Assume 4 stream for now

        !! Calculate Reflection coefficent for next layer
        !! Equation is: R(k+1,:,:) = T(k,:,:) * G(k,:,:) * R(k,:,:) * T(k,:,:) + R(k,:,:)
        R(k+1,:,:) = matmul(T(k,:,:),G(:,:))
        R(k+1,:,:) = matmul(R(k+1,:,:),R(k,:,:))
        R(k+1,:,:) = matmul(R(k+1,:,:),T(k,:,:))
        R(k+1,:,:) = R(k+1,:,:) + R(k,:,:)

        !! Calculate Transmission coefficent for next layer
        !! Equation is: T(k+1,:,:) = T(k,:,:) * G(k,:,:) * T(k,:,:)
        T(k+1,:,:) = matmul(T(k,:,:),G(:,:))
        T(k+1,:,:) = matmul(T(k+1,:,:),T(k,:,:))
      
        z(k+1,:) = (matmul(T(k,:,:),G(:,:)) - matmul(matmul(T(k,:,:),G(:,:)),R(k,:,:))) * (z(k,:) - gt*y(k,:)) + gt*y(k,:) + z(k,:)
        y(k+1,:) = (matmul(T(k,:,:),G(:,:)) + matmul(matmul(T(k,:,:),G(:,:)),R(k,:,:)) + I(:,:)) * y(k,:)

        sp(k+1,:) = y(k+1,:) * bm(k) + bl(k)*z(k+1,:)
        sm(k+1,:) = y(k+1,:) * bm(k) - bl(k)*z(k+1,:)

        cycle


      !else
        !! Must perform the doubling method for sub-layers until dtau_lim is reached
        !! Its probably best to split into equal dtau sub-layers

        !! Find number if doubling layers and initial sub-layer dtau
        nsub = int((log10(dtau(k)/(dtau_lim*uarr(1)))/log10(2.0_dp) + 1.0_dp))
        dtau_sub = dtau(k)/(2.0_dp**(nsub))

        gt = dtau_sub
  
        !! Do initial calculate at dtau_sub - take layer above as initial conditions

        !! Calculate Gamma 
        G_s(:,:) = I(:,:) - matmul(R(k,:,:),R(k,:,:))
        G_s(:,:) = matinv2(G_s(:,:)) ! Assume 4 stream for now

        !! Calculate Reflection coefficent for next layer
        !! Equation is: R(k+1,:,:) = T(k,:,:) * G(k,:,:) * R(k,:,:) * T(k,:,:) + R(k,:,:)
        R_s(:,:) = matmul(T(k,:,:),G(k,:,:))
        R_s(:,:) = matmul(R_s(:,:),R(k,:,:))
        R_s(:,:) = matmul(R_s(:,:),T(k,:,:))
        R_s(:,:) = R_s(:,:) + R(k,:,:)

        !! Calculate Transmission coefficent for next layer
        !! Equation is: T(k+1,:,:) = T(k,:,:) * G(k,:,:) * T(k,:,:)
        T_s(:,:) = matmul(T(k,:,:),G(k,:,:))
        T_s(:,:) = matmul(T_s(:,:),T(k,:,:))
      
        z_s(:) = (matmul(T(k,:,:),G(k,:,:)) - matmul(matmul(T(k,:,:),G(k,:,:)),R(k,:,:))) * (z(k,:) - gt*y(k,:)) + gt*y(k,:) + z(k,:)
        y_s(:) = (matmul(T(k,:,:),G(k,:,:)) + matmul(matmul(T(k,:,:),G(k,:,:)),R(k,:,:)) + I(:,:)) * y(k,:)

        sp_s(:) = y_s(:) * bm(k) + bl(k)*z_s(:)
        sm_s(:) = y_s(:) * bm(k) - bl(k)*z_s(:)

        !! We now have the values for the initial sub-layer

        do d = 2, nsub-1

          ! Increase sub layer dtau by 2
          dtau_sub = dtau_sub * 2.0_dp

          gt = dtau_sub

          ! Perform doubling


          ! Peform adding

        end do

        !! The values at the next computational level are the properties at the origonal level and last added sub-layer
        gt = dtau(k)

        !! Calculate Gamma 
        G(:,:) = I(:,:) - matmul(R_s(:,:),R_s(:,:))
        G(:,:) = matinv2(G(:,:)) ! Assume 4 stream for now

        !! Calculate Reflection coefficent for next layer
        !! Equation is: R(k+1,:,:) = T(k,:,:) * G(k,:,:) * R(k,:,:) * T(k,:,:) + R(k,:,:)
        R(k+1,:,:) = matmul(T_s(:,:),G(:,:))
        R(k+1,:,:) = matmul(R(k+1,:,:),R_s(:,:))
        R(k+1,:,:) = matmul(R(k+1,:,:),T_s(:,:))
        R(k+1,:,:) = R(k+1,:,:) + R_s(:,:)

        !! Calculate Transmission coefficent for next layer
        !! Equation is: T(k+1,:,:) = T(k,:,:) * G(k,:,:) * T(k,:,:)
        T(k+1,:,:) = matmul(T_s(:,:),G(:,:))
        T(k+1,:,:) = matmul(T(k+1,:,:),T_s(:,:))
      
        z(k+1,:) = (matmul(T_s(:,:),G(:,:)) - matmul(matmul(T_s(:,:),G(:,:)),R_s(:,:))) * (z_s(:) - gt_s*y_s(:)) + gt_s*y_s(:) + z_s(:)
        y(k+1,:) = (matmul(T_s(:,:),G(:,:)) + matmul(matmul(T_s(:,:),G(:,:)),R_s(:,:)) + I(:,:)) * y_s(:)

        sp(k+1,:) = y(k+1,:) * bm(k) + bl(k)*z(k+1,:)
        sm(k+1,:) = y(k+1,:) * bm(k) - bl(k)*z(k+1,:)


      !end if

    end do

    !! Now we know the Reflection and Transmission and source coefficents for each level
    !! Do the final adding step starting at the top of the atmosphere to get the source coefficent for the entire atmosphere
    do k = 2, nlev
      Rac(:,:) =  
      Rca(:,:) = 
    end do

    !! Downward and upward fluxes are calculated as quadrature
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp
    do m = 1, nmu
      flx_up(:) = flx_up(:) + sp_add(:,m)*w(m)
      flx_down(:) = flx_down(:) + sm_add(:,m)*w(m)
    end do


  end subroutine lw_adding_doubling

  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(dp), intent(in) :: A(2,2)   !! Matrix
    real(dp)             :: B(2,2)   !! Inverse matrix
    real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0_dp/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = detinv * A(1,1)

  end function matinv2

  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(dp), intent(in) :: A(4,4)   !! Matrix
    real(dp)             :: B(4,4)   !! Inverse matrix
    real(dp)             :: detinv, s0, s1, s2, s3, s4, s5, c5, c4, c3, c2, c1, c0

    s0 = A(1,1) * A(2,2) - A(2,1) * A(1,2)
    s1 = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    s2 = A(1,1) * A(2,4) - A(2,1) * A(1,4)
    s3 = A(1,2) * A(2,3) - A(2,2) * A(1,3)
    s4 = A(1,2) * A(2,4) - A(2,2) * A(1,4)
    s5 = A(1,3) * A(2,4) - A(2,3) * A(1,4)

    c5 = A(3,3) * A(4,4) - A(4,3) * A(3,4)
    c4 = A(3,2) * A(4,4) - A(4,2) * A(3,4)
    c3 = A(3,2) * A(4,3) - A(4,2) * A(3,3)
    c2 = A(3,1) * A(4,4) - A(4,1) * A(3,4)
    c1 = A(3,1) * A(4,3) - A(4,1) * A(3,3)
    c0 = A(3,1) * A(4,2) - A(4,1) * A(3,2)

    detinv = 1.0_dp / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)

    B(1,1) = ( A(2,2) * c5 - A(2,3) * c4 + A(2,4) * c3) * detinv
    B(1,2) = (-A(1,2) * c5 + A(1,3) * c4 - A(1,4) * c3) * detinv
    B(1,3) = ( A(4,2) * s5 - A(4,3) * s4 + A(4,4) * s3) * detinv
    B(1,4) = (-A(3,2) * s5 + A(3,3) * s4 - A(3,4) * s3) * detinv

    B(2,1) = (-A(2,1) * c5 + A(2,3) * c2 - A(2,4) * c1) * detinv
    B(2,2) = ( A(1,1) * c5 - A(1,3) * c2 + A(1,4) * c1) * detinv
    B(2,3) = (-A(4,1) * s5 + A(4,3) * s2 - A(4,4) * s1) * detinv
    B(2,4) = ( A(3,1) * s5 - A(3,3) * s2 + A(3,4) * s1) * detinv

    B(3,1) = ( A(2,1) * c4 - A(2,2) * c2 + A(2,4) * c0) * detinv
    B(3,2) = (-A(1,1) * c4 + A(1,2) * c2 - A(1,4) * c0) * detinv
    B(3,3) = ( A(4,1) * s4 - A(4,2) * s2 + A(4,4) * s0) * detinv
    B(3,4) = (-A(3,1) * s4 + A(3,2) * s2 - A(3,4) * s0) * detinv

    B(4,1) = (-A(2,1) * c3 + A(2,2) * c1 - A(2,3) * c0) * detinv
    B(4,2) = ( A(1,1) * c3 - A(1,2) * c1 + A(1,3) * c0) * detinv
    B(4,3) = (-A(4,1) * s3 + A(4,2) * s1 - A(4,3) * s0) * detinv
    B(4,4) = ( A(3,1) * s3 - A(3,2) * s1 + A(3,3) * s0) * detinv

  end function matinv4

  subroutine ludcmp(A,n,np,indx,D)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n), intent(out) :: indx
    real(dp), intent(out) :: D

    integer, parameter :: nmax = 100
    real(dp), parameter :: tiny = 1.0e-20_dp
    real(dp), dimension(nmax) :: vv

    integer :: i, j, k, imax
    real(dp) :: aamax, dum, sum

    D = 1.0_dp

    do i = 1, n
      aamax = 0.0_dp
      do j = 1, n
        if (abs(A(i,j)) > aamax) then
          aamax = abs(A(i,j))
        end if
      end do
      if (aamax == 0.0_dp) then
        print*, 'singualr matrix in LU decomp!'
        stop
      end if
      vv(i) = 1.0_dp/aamax
    end do

  
    do j = 1, n
      do i = 1, j-1
        sum = A(i,j)
        do k = 1, i-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
      end do
      aamax = 0.0_dp
      do i = j, n
        sum = A(i,j)
        do k = 1, j-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        end if
      end do
      if (j /= imax) then
        do k = 1, n
          dum = A(imax,k)
          A(imax,k) = A(j,k)
          A(j,k) = dum
        end do 
        D = -D
        vv(imax) = vv(j)
      end if
      indx(j) = imax
      if (A(j,j) <= tiny) then
        A(j,j) = tiny
      end if
      if (j /= n) then
        dum = 1.0_dp/A(j,j)
        do i = j+1, n
          A(i,j) = A(i,j)*dum
        end do
      end if
    end do

  end subroutine ludcmp

  subroutine lubksb(A, n, np, indx, B)
    implicit none

    integer, intent(in) :: n, np
    integer, dimension(n), intent(in) :: indx
    real(dp), dimension(np,np), intent(in) :: A

    real(dp), dimension(n), intent(out) :: B

    integer :: i, j, ii, ll
    real(dp) :: sum

    ii = 0

    do i = 1, n
      ll = indx(i)
      sum = B(ll)
      b(ll) = b(i)
      if (ii /= 0) then
        do j = ii,i-1
          sum = sum - A(i,j)*B(j)
        end do
      else if (sum /= 0.0_dp) then
        ii = i
      end if
      B(i) = sum
    end do

    do i = n, 1, -1
      sum = B(i)
      if (i < n) then
        do j = i+1, n
          sum = sum - A(i,j)*B(j)
        end do
      end if
      B(i) = sum/A(i,i)
    end do
    
  end subroutine lubksb

  subroutine inv_LU(A,n,np,Y)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n) :: indx
    real(dp), dimension(np,np), intent(out) :: Y

    real(dp) :: D
    integer :: i, j

    do i = 1, n
      do j = 1, n
        Y(i,j) = 0.0_dp
      end do
      Y(i,i) = 1.0_dp
    end do

    call ludcmp(A,n,np,indx,D)

    do j = 1, n
      call lubksb(A,n,np,indx,Y(1,j))
    end do

  end subroutine inv_LU

end module lw_AD_mod