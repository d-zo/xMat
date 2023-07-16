   ! --------------------------------------------------------------- !
   pure function Vec_To_Mat33(vec6) result(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(3, 3) :: mat33
      ! ------------------------------------------------------------ !
      mat33 = reshape([vec6(1), vec6(4), vec6(6), &
                       vec6(4), vec6(2), vec6(5), &
                       vec6(6), vec6(5), vec6(3)], [3, 3])
   end function Vec_To_Mat33
