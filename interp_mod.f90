module interp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  contains

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
  
    ! Perform linear interpolation in linear space
    subroutine linear_interp(xval, x1, x2, y1, y2, yval)
      implicit none
  
      real(dp), intent(in) :: xval, y1, y2, x1, x2
      real(dp), intent(out) :: yval
  
      yval = (y1 * (x2 - xval) + y2 * (xval - x1))/(x2 - x1)
  
    end subroutine linear_interp

    ! Perform bilinear interpolation in linear space
    subroutine bilinear_interp(xval, yval, x12, y12, matrix, P)
      implicit none

      real(dp), intent(in) :: xval, yval
      real(dp), dimension(2), intent(in) :: x12, y12
      real(dp), dimension(2,2), intent(in) :: matrix
      real(dp), intent(out) :: P
      real(dp) :: yval1, yval2

      yval1 = (matrix(1,1) * (x12(2) - xval) + matrix(2,1) * (xval - x12(1)))/(x12(2) - x12(1))
      yval2 = (matrix(1,2) * (x12(2) - xval) + matrix(2,2) * (xval - x12(1)))/(x12(2) - x12(1))

      P = ( 1.0_dp/((x12(2)-x12(1))*(y12(2)-y12(1))) ) * ( matrix(1,1)*(x12(2)-xval)*(y12(2)-yval) + &
                                                           matrix(2,1)*(xval-x12(1))*(y12(2)-yval) + &
                                                           matrix(1,2)*(x12(2)-xval)*(yval-y12(1)) + &
                                                           matrix(2,2)*(xval-x12(1))*(yval-y12(1)) )
    end subroutine bilinear_interp

end module interp_mod
