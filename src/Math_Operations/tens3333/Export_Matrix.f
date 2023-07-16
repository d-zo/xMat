   ! --------------------------------------------------------------- !
   pure function Export_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat_out may need to be adjusted hereafter
      real(dp), dimension(3, 3), intent(in) :: mat
      real(dp), dimension(num_dimensions+num_shear) :: mat_out
      ! ------------------------------------------------------------ !
      integer, dimension(4), parameter :: indexmap = [1, 2, 3, 1]
      integer :: idx

      mat_out = 0.0_dp
      do idx = 1, num_dimensions
         mat_out(idx) = mat(idx, idx)
      end do

      do idx = 1, num_shear
         mat_out(num_dimensions+idx) = 0.5_dp*(mat(indexmap(idx), indexmap(idx+1)) &
                                               + mat(indexmap(idx+1), indexmap(idx)))
      end do
   end function Export_Matrix
