   ! --------------------------------------------------------------- !
   pure function Norm(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp) :: Norm
      ! ------------------------------------------------------------ !
      Norm = sqrt(Double_Contraction22(vec9, vec9))                  ! Matrix 2-Norm (Frobenius norm) `||\mathbf{M}|| = \sqrt{\sum_i\sum_j|M_{ij}|^2}`
   end function Norm
