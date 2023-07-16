   ! --------------------------------------------------------------- !
   pure function Tensor_Partialtrace(mat66)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66
      real(dp) :: Tensor_Partialtrace
      ! ------------------------------------------------------------ !
      Tensor_Partialtrace = mat66(1, 1) + mat66(2, 2) + mat66(3, 3)
   end function Tensor_Partialtrace
