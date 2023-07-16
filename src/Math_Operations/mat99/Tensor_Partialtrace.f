   ! --------------------------------------------------------------- !
   pure function Tensor_Partialtrace(mat99)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: mat99
      real(dp) :: Tensor_Partialtrace
      ! ------------------------------------------------------------ !
      Tensor_Partialtrace = mat99(1, 1) + mat99(2, 2) + mat99(3, 3)
   end function Tensor_Partialtrace
