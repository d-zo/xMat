   ! --------------------------------------------------------------- !
   pure function Export_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat_out may need to be adjusted hereafter
      real(dp), dimension(6), intent(in) :: mat
      real(dp), dimension(num_dimensions+num_shear) :: mat_out
      ! ------------------------------------------------------------ !
      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      mat_out((num_dimensions + 1):(num_dimensions + num_shear)) = mat(4:(3 + num_shear))
   end function Export_Matrix
