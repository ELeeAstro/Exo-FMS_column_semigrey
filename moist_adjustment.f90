module moist_adj_mod
  use, intrinsic :: iso_fortran_env
  use phys, only: Rstar,H2O_MolecularWeight,H2O_cp
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  integer, parameter :: nb_convsteps = 10
  real(dp), parameter :: R = Rstar/H2O_MolecularWeight ! J/(kg*K) specific gas constant of water vapor
  real(dp), parameter :: Rcp = R/H2O_cp                ! cp in J/(kg*K) specific heat constant of water vapor
  contains

    !Dew point temperature
    function Tdew(p_dummy) result(res)
      use phys, only: H2O_L_vaporization,H2O_TriplePointT,H2O_TriplePointP,satvps
      implicit none
      real(dp), parameter :: L = H2O_L_vaporization        ! J/kg, latent heat of condensation of water vapor at 300K
                                                       ! should be H2O_L_sublimation if T < H2O_TriplePointT
      real(dp), parameter :: Tref = 350.0_dp           ! Reference temperature
      real(dp) :: pref                                 ! Saturation vapor pressure, arguments=T,T0,e0,MolecularWeight,LatentHeat

      real(dp), intent(in) :: p_dummy
      real(dp) :: res
      pref = satvps(Tref,H2O_TriplePointT,H2O_TriplePointP,H2O_MolecularWeight,L)
      res = Tref/(1.0_dp-(Tref*R/L)*log(p_dummy/pref))
    end function Tdew

    !Moist adjustment routine.
    subroutine moist_adj(Tmid,pmid,pint,t_step,dT_conv)
      implicit none
      real(dp), intent(in) :: t_step
      real(dp),dimension(:), intent(in) :: Tmid,pmid,pint ! pint isn't used yet, but could be adjusted here in case we
                                                               ! decide to account for surface condensation
      real(dp), dimension(:), intent(out) :: dT_conv
      real(dp), dimension(size(pmid)) :: Tmid_cc
      integer :: i, iter
      logical :: did_adj

      Tmid_cc = Tmid

      do iter = 1, nb_convsteps
        did_adj = .False.
        !------------------- L = cst -------------------
        do i=1,size(Tmid_cc) !Downward pass
          if (Tmid_cc(i) .lt. Tdew(pmid(i))) then
            Tmid_cc(i)=Tdew(pmid(i)) ! temperature stays the same during the phase change
            did_adj = .True.
          end if
        end do
        do i=size(Tmid_cc)-1,1,-1 !Upward pass
          if (Tmid_cc(i) .lt. Tdew(pmid(i))) then
            Tmid_cc(i)=Tdew(pmid(i))
            did_adj = .True.
          end if
        end do

        ! If no adjustment required, exit the loop
        if (did_adj .eqv. .False.) then
          exit
        end if

      end do

      ! Change in temperature is Tmid_cc - Tmid
      dT_conv(:) = (Tmid_cc(:) - Tmid(:))/t_step
      !print*, Tmid_cc(i), Tdew(pmid(i))

    end subroutine moist_adj

end module moist_adj_mod
