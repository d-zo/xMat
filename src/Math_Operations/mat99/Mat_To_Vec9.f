   ! --------------------------------------------------------------- !
   pure function Mat_To_Vec9(mat) result(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: mat
      real(dp), dimension(9) :: vec9
      ! ------------------------------------------------------------ !
      vec9 = mat
   end function Mat_To_Vec9
