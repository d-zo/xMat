   ! --------------------------------------------------------------- !
   pure function Squared(vec9) result(sq_vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(9) :: sq_vec9
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: temp_mat

      temp_mat = Vec_To_Mat33(vec9=vec9)
      temp_mat = matmul(temp_mat, temp_mat)                          ! Square `\mathbf{M}^2 = \mathbf{M}\cdot \mathbf{M}`
      sq_vec9 = Mat33_To_Vec(mat33=temp_mat)
   end function Squared
