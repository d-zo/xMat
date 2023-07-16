   ! --------------------------------------------------------------- !
   pure function Import_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat has to be adjusted beforehand
      real(dp), dimension(num_dimensions+num_shear), intent(in) :: mat
      real(dp), dimension(3, 3) :: mat_out
      ! ------------------------------------------------------------ !
      integer, dimension(4), parameter :: indexmap = [1, 2, 3, 1]
      integer :: idx

      mat_out = 0.0_dp
      do idx = 1, num_dimensions
         mat_out(idx, idx) = mat(idx)
      end do
      do idx = 1, num_shear
         mat_out(indexmap(idx), indexmap(idx+1)) = mat(num_dimensions + idx)
         mat_out(indexmap(idx+1), indexmap(idx)) = mat(num_dimensions + idx)
      end do
   end function Import_Matrix
