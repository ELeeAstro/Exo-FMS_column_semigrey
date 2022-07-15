module cloud_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64
  !! Dimensions of tables
  integer, parameter :: nb_r0 = 51
  integer, parameter :: nb_T  = 22

  real(dp), dimension(nb_r0,nb_T) :: aa_matrix, gg_matrix, kk_matrix
  real(dp), dimension(nb_r0)      :: radii
  real(dp), dimension(nb_T)       :: temperatures

  contains

    ! Computes the mean radius using the mixing ratio, density, and number density of the condensates
    subroutine mean_radius(q_c, N_c, sigma, r_0)
      use phys, only: pi, H2O_rho_liquid
      implicit none
      real(dp), intent(in) :: q_c, N_c, sigma
      real(dp), intent(out) :: r_0
      real(dp) :: rho_c

      rho_c = H2O_rho_liquid

      r_0 = ( (3.0_dp/(4.0_dp*pi))*(q_c/rho_c)*(1.0_dp/N_c)*exp(-((9.0_dp/2.0_dp)*sigma**2)) )**(1.0_dp/3.0_dp)

    end subroutine mean_radius

    ! Reads single scattering albedo, anisotropic factor, and opacity coefficient tables
    subroutine read_optical_properties()
      implicit none
      character(len=150) :: single_scattering_albedos, anisotropic_factors, opacity_coefficients, radius_temperature_table
      integer :: i
      
      single_scattering_albedos = './Cloud_parameters/H2O_rosselandMean_qscat.txt'
      anisotropic_factors       = './Cloud_parameters/H2O_rosselandMean_gg.txt'
      opacity_coefficients      = './Cloud_parameters/H2O_rosselandMean_qext.txt'
      radius_temperature_table  = './Cloud_parameters/H2O_rosselandMean_RTtable.txt'

      ! 51 radii and 22 temperatures
      ! Reading tables
      open(11,file=single_scattering_albedos)
      read(11,*)
      do i=1,nb_r0
        read(11,*) aa_matrix(i,:) ! aa_matrix(i,j), i is the radius index, j is the temperature index
      end do
      close(11)

      open(11,file=anisotropic_factors)
      read(11,*)
      do i=1,nb_r0
        read(11,*) gg_matrix(i,:) ! gg_matrix(i,j), i is the radius index, j is the temperature index
      end do
      close(11)

      open(11,file=opacity_coefficients)
      read(11,*)
      do i=1,nb_r0
        read(11,*) kk_matrix(i,:) ! kk_matrix(i,j), i is the radius index, j is the temperature index
      end do
      close(11)

      open(11,file=radius_temperature_table)
      read(11,*)
      read(11,*)
      read(11,*) radii
      read(11,*) temperatures
      close(11)

    end subroutine read_optical_properties

    ! Computes the single scattering albedo, anisotropic factor, and opacity coefficient through interpolation
    subroutine optical_properties(q_c, N_c, sigma, Tmid, aa, gg, kk)
      use interp_mod, only: bilinear_interp
      implicit none
      real(dp), intent(in) :: q_c, sigma, Tmid
      real(dp), intent(in) :: N_c
      real(dp), intent(out) :: aa, gg, kk
      real(dp) :: r_0, T_closest, T_closest_lower, T_closest_upper, r_closest, r_closest_lower, r_closest_upper
      character(len=150) :: single_scattering_albedos, anisotropic_factors, opacity_coefficients, radius_temperature_table
      real(dp), dimension(2) :: r2, T2
      real(dp), dimension(2,2) :: aa_matrix_22, gg_matrix_22, kk_matrix_22

      call mean_radius(q_c, N_c, sigma, r_0)

      r_closest = minloc(abs(radii-(radii/radii)*r_0), DIM=1) ! returns the index of the closest element to r_0(i) in radii
      if (radii(r_closest) < r_0) then ! lower neighbor
        r_closest_lower = r_closest
        r_closest_upper = r_closest+1
      else if (radii(r_closest) > r_0) then ! upper neighbor
        r_closest_upper = r_closest
        r_closest_lower = r_closest-1
      else if (radii(r_closest) == r_0) then
      end if

      T_closest = minloc(abs(temperatures-(temperatures/temperatures)*Tmid), DIM=1) 
      ! returns the index of the closest element to Tmid in temperatures

      if (temperatures(T_closest) < Tmid) then ! lower neighbor
        T_closest_lower = T_closest
        T_closest_upper = T_closest+1
      else if (temperatures(T_closest) > Tmid) then ! upper neighbor
        T_closest_upper = T_closest
        T_closest_lower = T_closest-1
      else if (temperatures(T_closest) == Tmid) then
      end if

      r2           = (/ radii(r_closest_lower), radii(r_closest_upper) /)
      T2           = (/ temperatures(T_closest_lower), temperatures(T_closest_upper) /)

      aa_matrix_22 = aa_matrix(r_closest_lower:r_closest_upper,T_closest_lower:T_closest_upper)
      call bilinear_interp(r_0,Tmid,r2,T2,aa_matrix_22,aa)

      gg_matrix_22 = gg_matrix(r_closest_lower:r_closest_upper,T_closest_lower:T_closest_upper)
      call bilinear_interp(r_0,Tmid,r2,T2,gg_matrix_22,gg)

      kk_matrix_22 = kk_matrix(r_closest_lower:r_closest_upper,T_closest_lower:T_closest_upper)
      call bilinear_interp(r_0,Tmid,r2,T2,kk_matrix_22,kk)

      if ( (Tmid .gt. maxval(temperatures)) .or. (Tmid .lt. minval(temperatures))&
      .or. (r_0 .gt. maxval(radii)) .or. (r_0 .lt. minval(radii)) ) then
        aa = aa_matrix(r_closest,T_closest) 
        gg = gg_matrix(r_closest,T_closest) 
        kk = kk_matrix(r_closest,T_closest) 
      end if

    end subroutine optical_properties

end module cloud_mod
