   ! --------------------------------------------------------------- !
   pure function Vec_To_Mat33(vec9) result(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(3, 3) :: mat33
      ! ------------------------------------------------------------ !
      mat33 = reshape([ &
         vec9(1), vec9(5), vec9(9), &
         vec9(4), vec9(2), vec9(7), &
         vec9(8), vec9(6), vec9(3)], [3, 3])
   end function Vec_To_Mat33
