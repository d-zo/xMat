   ! --------------------------------------------------------------- !
   pure function Dyadic_Product22(vec9a, vec9b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9a, vec9b
      real(dp), dimension(9, 9) :: Dyadic_Product22
      ! ------------------------------------------------------------ !
      Dyadic_Product22 = matmul(reshape(vec9a, [9, 1]), &            ! Dyadic product `T_{ijkl} = M_{A,ij}\cdot M_{B,kl}`
                                reshape(vec9b, [1, 9]))
   end function Dyadic_Product22
