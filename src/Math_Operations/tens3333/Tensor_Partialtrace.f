   ! --------------------------------------------------------------- !
   pure function Tensor_Partialtrace(tens3333)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tens3333
      real(dp) :: Tensor_Partialtrace
      ! ------------------------------------------------------------ !
      Tensor_Partialtrace = tens3333(1, 1, 1, 1) + tens3333(2, 2, 2, 2) + tens3333(3, 3, 3, 3)
   end function Tensor_Partialtrace
