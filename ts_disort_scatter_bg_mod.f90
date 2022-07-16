!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Dec 2021 : Bezier interpolation
! sw/lw: Two-stream DISORT version, modified by Xianyu Tan to include an internal heat source.
! Converted to band-grey by Ryan Boukrouche
!        Pros: Stable, performs accurate scattering calculations tried and tested, reliable model.
!        Cons: Slower than other methods.
!!!

module ts_disort_scatter_bg_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  public :: ts_disort_scatter_bg
  private :: linear_log_interp

contains

  subroutine ts_disort_scatter_bg(surf, Bezier, f_star, n_bands, wn_edges, nlay, nlev, Ts, Tl, &
             & pl, pe, tau_Ve, tau_IRe, tau_bg, mu_z, F0, Tint, AB, sw_a, sw_g, lw_a, &
             & lw_g, F_net, diffuse_up, diffuse_down, direct_beam, cff, scff, opr, asr)
    implicit none

    !! Additional input variables for band-grey
    integer, intent(in) :: n_bands
    real(dp), dimension(n_bands,2), intent(in) :: wn_edges

    !! Input variables
    logical, intent(in) :: surf, Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, mu_z, Tint, AB, Ts
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g
    real(dp), dimension(n_bands,nlev), intent(in) :: tau_bg
    real(dp), dimension(n_bands), intent(in) :: f_star

    !! Output variables
    real(dp), dimension(n_bands), intent(out) :: opr, asr
    real(dp), dimension(n_bands,nlev), intent(out) :: F_net, diffuse_up, diffuse_down, direct_beam
    real(dp), dimension(n_bands,nlay), intent(out) :: cff
    real(dp), dimension(n_bands), intent(out) :: scff

    !! Work variables
    real(dp), dimension(n_bands,nlev) :: lw_net, sw_net
    real(dp), dimension(n_bands,nlay) :: cff_lw, cff_sw
    real(dp), dimension(n_bands) :: scff_sw, opr_sw
    real(dp) :: umu0, fbeam
    logical :: planck
    integer :: i
    real(dp), dimension(nlev) :: Te

    !! Conversion arrays from FMS to DISORT dimensions
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(maxulv) :: pe_maxsize
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc, utau
    real(dp), dimension(maxcly) :: gg, ssalb
    real(dp), dimension(n_bands,maxulv) :: lw_net_maxsize, sw_net_maxsize, diffuse_lw_up_maxsize,&
                                           diffuse_sw_up_maxsize, diffuse_lw_down_maxsize, diffuse_sw_down_maxsize,&
                                           direct_beam_lw_maxsize, direct_beam_sw_maxsize 

    !! Additional variables for band-grey
    integer :: b

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

    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))

    if (surf .eqv. .True.) then
      Te(nlev) = Ts ! If surface, temperature at edge is the surface temperature
    else
      Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))
    end if

    Te_0(0:nlay) = Te(1:nlev)

    !! Shortwave flux calculation
    if (n_bands .eq. 1) then

      if (mu_z > 0.0_dp) then
        planck = .False.
        gg(1:nlay) = sw_g(:) ! gg, ssalb, dtauc, utau are size 200, net_F is size 201, diffuse_up(b,:) is size 55.
        ssalb(1:nlay) = sw_a(:)
        fbeam = (1.0_dp-AB) * F0
        umu0 = mu_z
        wvnmlo = 0.0_dp
        wvnmhi = 1.0e7_dp
        utau(1:nlev) = tau_Ve(:)
        do i = 1, nlay
          dtauc(i) = (tau_Ve(i+1) - tau_Ve(i))
        end do
        call CALL_TWOSTR (nlay,pe,Te_0,gg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,&
                          sw_net_maxsize(1,:),diffuse_sw_up_maxsize(1,:),diffuse_sw_down_maxsize(1,:),&
                          direct_beam_sw_maxsize(1,:),cff_sw(1,:),scff_sw(1),opr_sw(1))
        sw_net(1,:) = sw_net_maxsize(1,1:nlev)
        direct_beam(1,:) = direct_beam_sw_maxsize(1,1:nlev)
      else
        sw_net(1,:) = 0.0_dp
        direct_beam(1,:) = 0.0_dp
      end if

      !! Longwave two-stream flux calculation
      planck = .True.
      gg(1:nlay) = lw_g(:)
      ssalb(1:nlay) = lw_a(:)
      fbeam = 0.0_dp
      umu0 = 1.0_dp
      wvnmlo = 0.0_dp
      wvnmhi = 1.0e7_dp
      utau(1:nlev) = tau_IRe(:)
      do i = 1, nlay
        dtauc(i) = (tau_IRe(i+1) - tau_IRe(i))
      end do
      call CALL_TWOSTR (nlay,pe,Te_0,gg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,&
                        lw_net_maxsize(1,:),diffuse_lw_up_maxsize(1,:),diffuse_lw_down_maxsize(1,:),&
                        direct_beam_lw_maxsize(1,:),cff_lw(1,:),scff(1),opr(1))
                
      lw_net(1,:)  = lw_net_maxsize(1,1:nlev)

      !! Net fluxes at each level (includes the direct beam)
      F_net(1,:) = lw_net(1,1:nlev) + sw_net(1,1:nlev)

      !! Note sure how to calculate yet
      asr(1) = 0.0_dp

      !! Upward and downward fluxes
      diffuse_up(1,:)   = diffuse_lw_up_maxsize(1,1:nlev)   + diffuse_sw_up_maxsize(1,1:nlev)
      diffuse_down(1,:) = diffuse_lw_down_maxsize(1,1:nlev) + diffuse_sw_down_maxsize(1,1:nlev)
      direct_beam(1,:)  = direct_beam_lw_maxsize(1,1:nlev)  + direct_beam_sw_maxsize(1,1:nlev)

      !! CFF (already the right size)
      cff(1,:) = cff_lw(1,:) + cff_sw(1,:)

    else if (n_bands .ge. 4) then

      ! Direct beam calculation
      if (mu_z > 0.0_dp) then
        planck = .False.
        umu0 = mu_z
        do b = 1, n_bands
          gg(1:nlay) = sw_g(:)
          ssalb(1:nlay) = sw_a(:)
          fbeam = (1.0_dp-AB) * F0 * f_star(b)
          utau(1:nlev) = tau_bg(b,:)
          do i = 1, nlay
            dtauc(i) = tau_bg(b,i+1) - tau_bg(b,i)
          end do
          !print*, "pe = ", pe
          !print*, "tau_bg(b,:) = ", tau_bg(b,:)
          call CALL_TWOSTR (nlay,pe,Te_0,gg,ssalb,dtauc,nlev,utau,planck,wn_edges(b,2),wn_edges(b,1),Tint,&
                            fbeam,umu0,sw_net_maxsize(b,:),diffuse_sw_up_maxsize(b,:),diffuse_sw_down_maxsize(b,:),&
                            direct_beam_sw_maxsize(b,:),cff_sw(b,:),scff_sw(b),opr_sw(b))

          sw_net(b,:) = sw_net_maxsize(b,1:nlev)
          direct_beam(b,:) = direct_beam_sw_maxsize(b,1:nlev)

        end do
        !Subtract window fluxes from the IR band
        sw_net(1,:) = sw_net(1,:) - sw_net(2,:) - sw_net(3,:)
        direct_beam(1,:) = direct_beam(1,:) - direct_beam(2,:) - direct_beam(3,:)

      else 
        sw_net(:,:) = 0.0_dp
        direct_beam(:,:) = 0.0_dp
      end if
    
      ! Diffuse flux calculation
      planck = .True.
      umu0 = 1.0_dp
      do b = 1, n_bands
        gg(1:nlay) = lw_g(:)
        ssalb(1:nlay) = lw_a(:)
        fbeam = 0.0_dp
        utau(1:nlev) = tau_bg(b,:)
        do i = 1, nlay
          dtauc(i) = tau_bg(b,i+1) - tau_bg(b,i)
          !print*, "dtauc(i), tau_bg(b,i+1), tau_bg(b,i) = ", dtauc(i), tau_bg(b,i+1), tau_bg(b,i)
        end do
        call CALL_TWOSTR (nlay,pe,Te_0,gg,ssalb,dtauc,nlev,utau,planck,wn_edges(b,2),wn_edges(b,1),Tint,&
                          fbeam,umu0,lw_net_maxsize(b,:),diffuse_lw_up_maxsize(b,:),diffuse_lw_down_maxsize(b,:),&
                          direct_beam_lw_maxsize(b,:),cff_lw(b,:),scff(b),opr(b))
                    
        !print*, "cff_lw(b,:) = ", cff_lw(b,:)
        !! Net fluxes at each level
        lw_net(b,:)  = lw_net_maxsize(b,1:nlev)

        !! Net fluxes at each level (includes the direct beam)
        F_net(b,:) = lw_net(b,1:nlev) + sw_net(b,1:nlev)

        !! Note sure how to calculate yet
        asr(b) = 0.0_dp

        !! Upward and downward fluxes
        diffuse_up(b,:)   = diffuse_lw_up_maxsize(b,1:nlev)   + diffuse_sw_up_maxsize(b,1:nlev)
        diffuse_down(b,:) = diffuse_lw_down_maxsize(b,1:nlev) + diffuse_sw_down_maxsize(b,1:nlev)
        direct_beam(b,:)  = direct_beam_lw_maxsize(b,1:nlev)  + direct_beam_sw_maxsize(b,1:nlev)

        print*, "Fup",b,"1), sub, wn_edges(",b,",2),wn_edges(",b,",1) = ", diffuse_up(b,1),&
        diffuse_up(1,1)-diffuse_up(2,1)-diffuse_up(3,1), wn_edges(b,2),wn_edges(b,1)

        !! CFF (already the right size)
        cff(b,:) = cff_lw(b,:) + cff_sw(b,:)

        ! Checking for NaNs
        do i = 1, nlay+1
          if (isnan(pe(i))) stop '"pe(i)" is a NaN'
          if (isnan(Te(i))) stop '"Te(i)" is a NaN'
          
        end do
        do i = 1, nlay
          if (isnan(gg(i))) stop '"gg(i)" is a NaN'
          if (isnan(ssalb(i))) stop '"ssalb(i)" is a NaN'
          if (isnan(dtauc(i))) stop '"dtauc(i)" is a NaN'
        end do
        do i = 1, nlay+1
          if (isnan(utau(i))) stop '"utau(i)" is a NaN'
        end do
        if (isnan(wn_edges(b,2))) stop '"wn_edges(b,2)" is a NaN'
        if (isnan(wn_edges(b,1))) stop '"wn_edges(b,1)" is a NaN'
        if (isnan(Tint)) stop '"Tint" is a NaN'
        if (isnan(fbeam)) stop '"fbeam" is a NaN'
        if (isnan(umu0)) stop '"umu0" is a NaN'
        do i = 1, nlay+1
          if (isnan(F_net(b,i))) stop '"diffuse_net(b,i)" is a NaN'
          if (isnan(diffuse_up(b,i))) stop '"diffuse_up(b,i)" is a NaN'
          if (isnan(diffuse_down(b,i))) stop '"diffuse_up(b,i)" is a NaN'
          if (isnan(direct_beam(b,i))) stop '"direct_beam(b,i)" is a NaN'
        end do
        do i = 1, nlay
          if (isnan(cff(b,i))) stop '"cff(b,i)" is a NaN'
        end do
        if (isnan(scff(b))) stop '"scff(b)" is a NaN'
        if (isnan(opr(b))) stop '"opr(b)" is a NaN'
        if (isnan(asr(b))) stop '"asr(b)" is a NaN'

      end do

      !Subtract window fluxes from the IR band
      lw_net(1,:)       = lw_net(1,:)       - lw_net(2,:)       - lw_net(3,:)
      sw_net(1,:)       = sw_net(1,:)       - sw_net(2,:)       - sw_net(3,:)
      !print*, "diffuse_up(1,1), diffuse_up(2,1), diffuse_up(3,1), sub = ", diffuse_up(1,1), diffuse_up(2,1), diffuse_up(3,1),&
      diffuse_up(1,1)   - diffuse_up(2,1)   - diffuse_up(3,1)
      diffuse_up(1,:)   = diffuse_up(1,:)   - diffuse_up(2,:)   - diffuse_up(3,:)
      diffuse_down(1,:) = diffuse_down(1,:) - diffuse_down(2,:) - diffuse_down(3,:)
      cff(1,:)          = cff(1,:)          - cff(2,:)          - cff(3,:)
      scff(1)           = scff(1)           - scff(2)           - scff(3)
      opr(1)            = opr(1)            - opr(2)            - opr(3)
      asr(1)            = asr(1)            - asr(2)            - asr(3)

    end if

  end subroutine ts_disort_scatter_bg

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

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

  end subroutine bezier_interp

end module ts_disort_scatter_bg_mod
