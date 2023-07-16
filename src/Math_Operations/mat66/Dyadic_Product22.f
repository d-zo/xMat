   ! --------------------------------------------------------------- !
   pure function Dyadic_Product22(vec6a, vec6b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6a, vec6b
      real(dp), dimension(6, 6) :: Dyadic_Product22
      ! ------------------------------------------------------------ !
      Dyadic_Product22 = matmul(reshape(vec6a, [6, 1]), &            ! Dyadic product `T_{ijkl} = M_{A,ij}\cdot M_{B,kl}`
                                reshape(vec6b, [1, 6]))
   end function Dyadic_Product22
