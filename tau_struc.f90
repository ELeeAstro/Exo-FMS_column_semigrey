module tau_struc_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64
  contains
    
    subroutine tau_struct(nlay,grav,pint,kRoss,tau_edges)
      implicit none

      integer, intent(in) :: nlay
      real(dp), intent(in) :: grav
      real(dp), dimension(nlay+1), intent(in) :: pint
      real(dp), dimension(nlay), intent(in) :: kRoss
      real(dp), dimension(nlay+1), intent(out) :: tau_edges

      real(dp) :: tau_sum, tau_lay, deltaP, pref_broad
      integer :: k

      ! Reference pressure for pressure-broadened absorption, set at 0.1 bar
      pref_broad = 1e4_dp

      ! Initializing the running sum of optical depth
      tau_sum = 0.0_dp
      tau_edges(1) = tau_sum ! top edge of the atmosphere has zero optical depth

      ! Integrate from top to bottom
      do k = 1, nlay
        ! Pressure difference between layer edges
        deltaP = (pint(k+1) - pint(k))

        ! Optical depth of layer, assuming hydrostatic equilibrium
        tau_lay = (kRoss(k) * deltaP) / grav
        !print*, kRoss(k)

        ! Add to running sum
        tau_sum = tau_sum + tau_lay

        ! Optical depth structure is running sum
        tau_edges(k+1) = tau_sum

      end do

    end subroutine tau_struct

end module tau_struc_mod
