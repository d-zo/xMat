   ! --------------------------------------------------------------- !
   pure function Matrix_Rotation(vec6, rot_mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(3, 3), intent(in) :: rot_mat33
      real(dp), dimension(6) :: Matrix_Rotation
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: M_mat

      M_mat = Vec_To_Mat33(vec6=vec6)
      M_mat = matmul(rot_mat33, matmul(M_mat, transpose(rot_mat33))) ! Rotation `\mathbf{D}\cdot \mathbf{X}\cdot \mathbf{D}^T` of symmetric matrix `\mathbf{X}`
      Matrix_Rotation = Mat33_To_Vec(mat33=M_mat)                    ! results in another symmetric matrix
   end function Matrix_Rotation
