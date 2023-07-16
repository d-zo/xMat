   ! --------------------------------------------------------------- !
   pure function Vec9_To_Mat(vec9) result(mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(6) :: mat
      ! ------------------------------------------------------------ !
      mat = [vec9(1:3), 0.5_dp*(vec9(4)+vec9(5)), 0.5_dp*(vec9(6)+vec9(7)), 0.5_dp*(vec9(8)+vec9(9))]
   end function Vec9_To_Mat
