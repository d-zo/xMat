   ! --------------------------------------------------------------- !
   pure function Double_Contraction42(mat66, vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: Double_Contraction42
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: temp_vec

      temp_vec(1:3) = vec6(1:3)
      temp_vec(4:6) = 2.0_dp*vec6(4:6)
      Double_Contraction42 = matmul(mat66, temp_vec)                 ! Double contraction `\mathcal{T} : \mathbf{M} = \sum_{kl} T_{ijkl}\cdot M_{kl}`
   end function Double_Contraction42
