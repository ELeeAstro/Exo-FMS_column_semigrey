!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!
! sw: 
!     Pros: 
!     Cons: 
!!!


module sw_Feautrier_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_Feautrier_method
  public :: sw_Feautrier

contains

  subroutine sw_Feautrier(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev) :: mu_z
    real(dp) :: Finc, a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g


    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, nb
        sw_down_b(b,:) = 0.0_dp
        sw_up_b(b,:) = 0.0_dp
        if (Finc(b) <= 0.0_dp) then
          cycle
        end if
        do g = 1, ng
          call sw_Feautrier_method(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), & 
            & sw_down_g(g,:), sw_up_g(g,:))
          sw_down_b(b,:) = sw_down_b(b,:) + sw_down_g(g,:) * gw(g)
          sw_up_b(b,:) = sw_up_b(b,:) + sw_up_g(g,:) * gw(g)
        end do
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_Feautrier


  subroutine sw_Feautrier_method(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    integer :: k
    real(dp), dimension(nlev) :: cum_trans

    integer, parameter :: na = 1
    real(dp), dimension(nlev,na,na) :: A, B, C
    real(dp), dimension(nlev,na) :: L
    real(dp), dimension(nlev, na) :: j


    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        do k = 1, nlev
          flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
        end do
      end if

      flx_down(nlev) = flx_down(nlev) * (1.0_dp - w_surf_in) ! The surface flux for surface heating is the amount of flux absorbed by surface
      flx_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo

      return

    end if

    A(:,:,:) = 0.0_dp
    B(:,:,:) = 0.0_dp
    C(:,:,:) = 0.0_dp

    do k = 1, nlev

      A(k,1,1) = 

 
    end do

  end subroutine sw_Feautrier_method

end module sw_Feautrier_mod
