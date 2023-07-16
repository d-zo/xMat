   ! --------------------------------------------------------------- !
   pure function Double_Contraction22(vec9a, vec9b) result(vec_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9a, vec9b
      real(dp) :: vec_out
      ! ------------------------------------------------------------ !
      vec_out = dot_product(vec9a, vec9b)                            ! Double Contraction `\mathbf{M}_{A} : \mathbf{M}_{B} = \sum_{ij} M_{A,ij}\cdot M_{B,ij}`
   end function Double_Contraction22
