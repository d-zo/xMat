   ! --------------------------------------------------------------- !
   pure function Matrix_Rotation(mat, rot_mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat
      real(dp), dimension(3, 3), intent(in) :: rot_mat33
      real(dp), dimension(3, 3) :: Matrix_Rotation
      ! ------------------------------------------------------------ !
      Matrix_Rotation = matmul(rot_mat33, matmul(mat, transpose(rot_mat33)))
   end function Matrix_Rotation
