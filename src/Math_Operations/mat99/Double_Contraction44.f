   ! --------------------------------------------------------------- !
   pure function Double_Contraction44(mat99a, mat99b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: mat99a, mat99b
      real(dp), dimension(9, 9) :: Double_Contraction44
      ! ------------------------------------------------------------ !
      Double_Contraction44 = matmul(mat99a, mat99b)                  ! Double contraction `\mathcal{T} : \mathcal{F} = \sum_{kl} T_{ijkl}\cdot F_{klmn}`
   end function Double_Contraction44
