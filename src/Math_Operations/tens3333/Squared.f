   ! --------------------------------------------------------------- !
   pure function Squared(mat33) result(sq_mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp), dimension(3, 3) :: sq_mat33
      ! ------------------------------------------------------------ !
      sq_mat33 = matmul(mat33, mat33)                                ! Square `\mathbf{M}^2 = \mathbf{M}\cdot \mathbf{M}`
   end function Squared
