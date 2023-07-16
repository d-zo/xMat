   ! --------------------------------------------------------------- !
   pure function Import_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat has to be adjusted beforehand
      real(dp), dimension(num_dimensions+num_shear), intent(in) :: mat
      real(dp), dimension(9) :: mat_out
      ! ------------------------------------------------------------ !
      integer :: idx

      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      do idx = 1, num_shear
         mat_out(2 + 2*idx) = mat(num_dimensions + idx)
         mat_out(3 + 2*idx) = mat(num_dimensions + idx)
      end do
   end function Import_Matrix
