   ! --------------------------------------------------------------- !
   pure function Double_Contraction22(vec6a, vec6b) result(vec_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6a, vec6b
      real(dp) :: vec_out
      ! ------------------------------------------------------------ !
      vec_out = dot_product(vec6a, [vec6b(1:3), 2.0_dp*vec6b(4:6)])  ! Double Contraction `\mathbf{M}_{A} : \mathbf{M}_{B} = \sum_{ij} M_{A,ij}\cdot M_{B,ij}`
   end function Double_Contraction22
