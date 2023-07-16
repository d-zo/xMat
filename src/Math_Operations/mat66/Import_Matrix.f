   ! --------------------------------------------------------------- !
   pure function Import_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat has to be adjusted beforehand
      real(dp), dimension(num_dimensions+num_shear), intent(in) :: mat
      real(dp), dimension(6) :: mat_out
      ! ------------------------------------------------------------ !
      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      mat_out(4:(3 + num_shear)) = mat((num_dimensions + 1):(num_dimensions + num_shear))
   end function Import_Matrix
