module read_planck_mod 
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64
  real(dp), dimension(:,:), allocatable :: Planck_table

  contains

    subroutine read_planck(n_bands)
      implicit none
      integer, intent(in) :: n_bands
      character (len=20) :: n_bands_string
      integer :: i
      character(len=150) :: path_to_tables
      
      write(n_bands_string, '(I0.1)') n_bands

      path_to_tables = './planck_tables'
      ! 1st column is temperature, 1st row is headings
      open(11,file=trim(path_to_tables)//'/planck_table_'//trim(n_bands_string)//'bands_TPLKAVG.dat') 
      read(11,*) ! skip 1st row
      do i = 1, size(Planck_table,2)
        if (n_bands .eq. 1) then
          read(11,*) Planck_table(1,i), Planck_table(2,i)
        elseif (n_bands .eq. 4) then
          read(11,*) Planck_table(1,i), Planck_table(2,i), Planck_table(3,i), Planck_table(4,i), Planck_table(5,i)
        elseif (n_bands .eq. 5) then
          read(11,*) Planck_table(1,i), Planck_table(2,i), Planck_table(3,i), Planck_table(4,i), Planck_table(5,i),&
                     Planck_table(6,i)
        elseif (n_bands .eq. 6) then
          read(11,*) Planck_table(1,i), Planck_table(2,i), Planck_table(3,i), Planck_table(4,i), Planck_table(5,i),&
                     Planck_table(6,i), Planck_table(7,i)
        end if
      end do
      close(11)

    end subroutine read_planck

end module read_planck_mod
