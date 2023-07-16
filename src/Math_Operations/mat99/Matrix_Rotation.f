   ! --------------------------------------------------------------- !
   pure function Matrix_Rotation(vec9, rot_mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(3, 3), intent(in) :: rot_mat33
      real(dp), dimension(9) :: Matrix_Rotation
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: M_mat

      M_mat = Vec_To_Mat33(vec9=vec9)
      M_mat = matmul(rot_mat33, matmul(M_mat, transpose(rot_mat33)))
      Matrix_Rotation = Mat33_To_Vec(mat33=M_mat)
   end function Matrix_Rotation
