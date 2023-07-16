   ! --------------------------------------------------------------- !
   pure function Norm(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp) :: Norm
      ! ------------------------------------------------------------ !
      Norm = sqrt(Double_Contraction22(mat33, mat33))                ! Matrix 2-Norm (Frobenius norm) `||\mathbf{M}|| = \sqrt{\sum_i\sum_j|M_{ij}|^2}`
   end function Norm
