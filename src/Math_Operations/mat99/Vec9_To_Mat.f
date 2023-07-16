   ! --------------------------------------------------------------- !
   pure function Vec9_To_Mat(vec9) result(mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(9) :: mat
      ! ------------------------------------------------------------ !
      mat = vec9
   end function Vec9_To_Mat
