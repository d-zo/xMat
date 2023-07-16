   ! --------------------------------------------------------------- !
   pure function Determinant(vec6) result(det6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: det6
      ! ------------------------------------------------------------ !
      if (global_num_direct_components == 3) then
         det6 = vec6(1)*(vec6(2)*vec6(3) - vec6(5)*vec6(5)) &
              + vec6(4)*(vec6(5)*vec6(6) - vec6(3)*vec6(4)) &
              + vec6(6)*(vec6(4)*vec6(5) - vec6(2)*vec6(6))
      else
         det6 = vec6(1)*vec6(2) - vec6(4)*vec6(4)
      end if
   end function Determinant
