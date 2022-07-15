module spectral_partitioner_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64
  real(dp), dimension(:,:), allocatable :: spectrum

  contains

    subroutine read_stellar_spectrum(star)
      ! Reads the spectrum of a given star and converts wavelengths to m where needed.
      implicit none
      character(len=*), intent(in) :: star
      character(len=150) :: spectrum_path
      integer :: i
      select case(star)
        case('Sun')
          allocate(spectrum(2,3780))      
          spectrum_path = './Sun/lean_12' 
          open(11,file=spectrum_path)
          read(11,*)
          read(11,*)
          read(11,*)
          read(11,*)
          read(11,*)
          do i=1, size(spectrum,2)
            read(11,*) spectrum(1,i), spectrum(2,i)
          end do
          close(11) 

        case('GJ551')
          allocate(spectrum(2,54997))
          spectrum_path = './GJ551/sflux-GJ551.txt'
          open(11,file=spectrum_path)
          read(11,*)
          do i=1, size(spectrum,2)
            read(11,*) spectrum(1,i), spectrum(2,i) 
          end do
          close(11)
          spectrum(1,:) = spectrum(1,:) * 1e-9
          spectrum(2,:) = spectrum(2,:) * 1e3 ! [Y erg/cm^2/s/Hz] = 1000 * [X W/m^2/Hz]
        end select
    end subroutine read_stellar_spectrum

    subroutine spectral_partition(n_bands,star,f_star,I_star) 
      ! Given a number of spectral bands n_bands, returns the fraction of incident stellar flux in each band.
      implicit none
      integer, intent(in) :: n_bands
      character(len=*), intent(in) :: star
      real(dp), dimension(n_bands), intent(out) :: f_star, I_star
      real(dp), dimension(2) :: IR1, IR2, IR3, W1, W2, SW, UV, VIS, VIS1, VIS2
      real(dp) :: sum_IR1, sum_IR2, sum_IR3, sum_W1, sum_W2, sum_SW, sum_UV, sum_VIS, sum_VIS1, sum_VIS2
      integer :: i

      IR1  = (/ 800.00d-9, 3448.28d-9 /)
      IR2  = (/ 4545.45d-9, 7692.31d-9 /) 
      IR3  = (/ 20000.00d-9, 10000000.00d-9 /)
      W1   = (/ 3448.28d-9, 4545.45d-9 /)
      W2   = (/ 7692.31d-9, 20000.00d-9 /)
      SW   = (/ 246.91d-9, 800.00d-9 /)   
      UV   = (/ 246.91d-9, 400.00d-9 /)
      VIS  = (/ 400.00d-9, 800.00d-9 /)
      VIS1 = (/ 400.00d-9, 533.33d-9 /)
      VIS2 = (/ 533.33d-9, 800.00d-9 /)

      call read_stellar_spectrum(star)

      sum_IR1 = 0.0 ; sum_IR2 = 0.0 ; sum_IR3 = 0.0 ; sum_W1 = 0.0 ; sum_W2 = 0.0 
      sum_SW = 0.0 ; sum_UV = 0.0 ; sum_VIS = 0.0 ; sum_VIS1 = 0.0 ; sum_VIS2 = 0.0
      do i = 1, size(spectrum,2)
        if ((spectrum(1,i) .ge. IR1(1)) .and. (spectrum(1,i) .le. IR1(2))) then
          sum_IR1 = sum_IR1 + spectrum(2,i)
        end if
        if ((spectrum(1,i) .ge. IR2(1)) .and. (spectrum(1,i) .le. IR2(2))) then
          sum_IR2 = sum_IR2 + spectrum(2,i)
        end if
        if ((spectrum(1,i) .ge. IR3(1)) .and. (spectrum(1,i) .le. IR3(2))) then
          sum_IR3 = sum_IR3 + spectrum(2,i)
        end if
        I_star(1) = sum_IR1+sum_IR2+sum_IR3
        f_star(1) = I_star(1)/sum(spectrum(2,:)) 
  
        if ((spectrum(1,i) .ge. W1(1)) .and. (spectrum(1,i) .le. W1(2))) then
          sum_W1 = sum_W1 + spectrum(2,i)
        end if
        I_star(2) = sum_W1
        f_star(2) = I_star(2)/sum(spectrum(2,:))
  
        if ((spectrum(1,i) .ge. W2(1)) .and. (spectrum(1,i) .le. W2(2))) then
          sum_W2 = sum_W2 + spectrum(2,i)
        end if
        I_star(3) = sum_W2
        f_star(3) = I_star(3)/sum(spectrum(2,:))
  
        if (n_bands .eq. 4) then
  
          if ((spectrum(1,i) .ge. SW(1)) .and. (spectrum(1,i) .le. SW(2))) then
            sum_SW = sum_SW + spectrum(2,i)
          end if
          I_star(4) = sum_SW
          f_star(4) = I_star(4)/sum(spectrum(2,:))
  
        elseif (n_bands .eq. 5) then
  
          if ((spectrum(1,i) .ge. UV(1)) .and. (spectrum(1,i) .le. UV(2))) then
            sum_UV = sum_UV + spectrum(2,i)
          end if
          I_star(4) = sum_UV
          f_star(4) = I_star(4)/sum(spectrum(2,:))
  
          if ((spectrum(1,i) .ge. VIS(1)) .and. (spectrum(1,i) .le. VIS(2))) then
            sum_VIS = sum_VIS + spectrum(2,i)
          end if
          I_star(5) = sum_VIS
          f_star(5) = I_star(5)/sum(spectrum(2,:))
  
        elseif (n_bands .eq. 6) then
  
          if ((spectrum(1,i) .ge. UV(1)) .and. (spectrum(1,i) .le. UV(2))) then
            sum_UV = sum_UV + spectrum(2,i)
          end if
          I_star(4) = sum_UV
          f_star(4) = I_star(4)/sum(spectrum(2,:))
  
          if ((spectrum(1,i) .ge. VIS1(1)) .and. (spectrum(1,i) .le. VIS1(2))) then
            sum_VIS1 = sum_VIS1 + spectrum(2,i)
          end if
          I_star(5) = sum_VIS1
          f_star(5) = I_star(5)/sum(spectrum(2,:))
  
          if ((spectrum(1,i) .ge. VIS2(1)) .and. (spectrum(1,i) .le. VIS2(2))) then
            sum_VIS2 = sum_VIS2 + spectrum(2,i)
          end if
          I_star(6) = sum_VIS2
          f_star(6) = I_star(6)/sum(spectrum(2,:))
  
        end if
      end do

      deallocate(spectrum)

    end subroutine spectral_partition

end module spectral_partitioner_mod
