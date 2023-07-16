   ! --------------------------------------------------------------- !
   pure function Double_Contraction42(mat99, vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: mat99
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(9) :: Double_Contraction42
      ! ------------------------------------------------------------ !
      Double_Contraction42 = matmul(mat99, vec9)                     ! Double contraction `\mathcal{T} : \mathbf{M} = \sum_{kl} T_{ijkl}\cdot M_{kl}`
   end function Double_Contraction42
