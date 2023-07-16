   ! --------------------------------------------------------------- !
   pure function Double_Contraction22(mat33a, mat33b) result(res)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33a, mat33b
      real(dp) :: res
      ! ------------------------------------------------------------ !
      res = dot_product(reshape(mat33a, [9]), reshape(mat33b, [9]))  ! Double Contraction `\mathbf{M}_{A} : \mathbf{M}_{B} = \sum_{ij} M_{A,ij}\cdot M_{B,ij}`
   end function Double_Contraction22
