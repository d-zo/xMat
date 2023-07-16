   ! --------------------------------------------------------------- !
   pure function Mat_To_Vec9(mat) result(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: mat
      real(dp), dimension(9) :: vec9
      ! ------------------------------------------------------------ !
      vec9 = [mat(1:4), mat(4), mat(5), mat(5), mat(6), mat(6)]
   end function Mat_To_Vec9
