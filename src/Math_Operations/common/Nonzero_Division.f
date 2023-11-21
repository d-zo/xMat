   ! --------------------------------------------------------------- !
   elemental function Nonzero_Division(val, fac)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: val
      real(dp), intent(in) :: fac
      real(dp) :: Nonzero_Division
      ! ------------------------------------------------------------ !
      if (abs(fac) < setting_epsilon_extra) then
         Nonzero_Division = val
      else
         Nonzero_Division = val/fac
      end if
   end function Nonzero_Division
