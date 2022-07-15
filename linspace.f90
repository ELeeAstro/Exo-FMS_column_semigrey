module linspace_module
  implicit none
  contains
  ! Creates a 1-dimensional array with evenly spaced elements

    subroutine linspace(x,x_start,x_end,x_len)
      implicit none
      real(kind=8), dimension(:), intent(out) :: x
      real(kind=8), intent(in) :: x_start,x_end
      integer, intent(in) :: x_len
      real(kind=8) :: dx
      integer :: i
      dx=(x_end-x_start)/(x_len-1)
      x(1:x_len)=[(x_start + ((i-1)*dx), i=1,x_len)]
    end subroutine linspace

    function arange(x_start, x_end, increment) RESULT(res)
    ! Returns an array of reals given x_start,  x_end,  and increment values.
    ! Increment defaults to 1 if not provided.
      implicit none
      real(kind=8), intent(in) :: x_start,x_end
      real(kind=8), intent(in), optional :: increment
      real(kind=8), dimension(:), allocatable :: res
      real(kind=8) :: incr
      integer :: i, length
      if(present(increment)) then
        incr = increment
      else
        incr = 1
      end if
      length = (x_end-x_start+0.5D0*incr)/incr+1
      allocate(res(length))
      do concurrent(i = 1:length)
        res(i) = x_start+(i-1)*incr
      end do
    end function arange

END MODULE linspace_module
