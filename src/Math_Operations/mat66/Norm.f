   ! --------------------------------------------------------------- !
   pure function Norm(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: Norm
      ! ------------------------------------------------------------ !
      Norm = sqrt(Double_Contraction22(vec6, vec6))                  ! Matrix 2-Norm (Frobenius norm) `||\mathbf{M}|| = \sqrt{\sum_i\sum_j|M_{ij}|^2}`
   end function Norm
