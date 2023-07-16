   ! --------------------------------------------------------------- !
   pure function Value_In_Interval(value, limits)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: value
      real(dp), dimension(2), intent(in) :: limits
      logical :: Value_In_Interval
      ! ------------------------------------------------------------ !
      Value_In_Interval = .False.
      if ((value >= minval(limits)) .and. (value <= maxval(limits))) then
         Value_In_Interval = .True.
      end if
   end function Value_In_Interval
