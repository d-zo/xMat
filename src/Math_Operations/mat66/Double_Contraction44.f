   ! --------------------------------------------------------------- !
   pure function Double_Contraction44(mat66a, mat66b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66a, mat66b
      real(dp), dimension(6, 6) :: Double_Contraction44
      ! ------------------------------------------------------------ !
      real(dp), dimension(6, 6) :: temp_mat66a

      temp_mat66a(:, 1:3) = mat66a(:, 1:3)
      temp_mat66a(:, 4:6) = 2.0_dp*mat66a(:, 4:6)
      Double_Contraction44 = matmul(temp_mat66a, mat66b)             ! Double contraction `\mathcal{T} : \mathcal{F} = \sum_{kl} T_{ijkl}\cdot F_{klmn}`
   end function Double_Contraction44
