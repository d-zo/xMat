   ! --------------------------------------------------------------- !
   pure function Export_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat_out may need to be adjusted hereafter
      real(dp), dimension(9), intent(in) :: mat
      real(dp), dimension(num_dimensions+num_shear) :: mat_out
      ! ------------------------------------------------------------ !
      integer :: idx

      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      do idx = 1, num_shear
         mat_out(num_dimensions+idx) = 0.5_dp*(mat(2 + 2*idx)+mat(3 + 2*idx))
      end do
   end function Export_Matrix
