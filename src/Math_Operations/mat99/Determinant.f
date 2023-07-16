   ! --------------------------------------------------------------- !
   pure function Determinant(vec9) result(det9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp) :: det9
      ! ------------------------------------------------------------ !
      if (global_num_direct_components == 3) then
         det9 = vec9(1)*(vec9(2)*vec9(3) - vec9(6)*vec9(7)) &
              + vec9(4)*(vec9(6)*vec9(9) - vec9(3)*vec9(5)) &
              + vec9(8)*(vec9(5)*vec9(7) - vec9(2)*vec9(9))
      else
         det9 = vec9(1)*vec9(2) - vec9(4)*vec9(5)
      end if
   end function Determinant
