   ! --------------------------------------------------------------- !
   elemental function Sqrt_Of_Sum_Of_Squares(a, b)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: a, b
      real(dp) :: Sqrt_Of_Sum_Of_Squares
      ! ------------------------------------------------------------ !
      real(dp) :: abs_value_a, abs_value_b

      abs_value_a = abs(a)
      abs_value_b = abs(b)

      if (abs_value_b < setting_epsilon_extra) then
         Sqrt_Of_Sum_Of_Squares = 0.0_dp
      else if (abs_value_a > abs_value_b) then
         Sqrt_Of_Sum_Of_Squares = abs_value_a*sqrt(1.0_dp + (abs_value_b/abs_value_a)**2)
      else
         Sqrt_Of_Sum_Of_Squares = abs_value_b*sqrt(1.0_dp + (abs_value_a/abs_value_b)**2)
      end if
   end function Sqrt_Of_Sum_Of_Squares
