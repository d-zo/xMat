   ! --------------------------------------------------------------- !
   pure function Is_First_Call(step, iteration)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: step, iteration
      logical :: Is_First_Call
      ! ------------------------------------------------------------ !
      if ((step == 1) .and. (iteration == 1)) then
         Is_First_Call = .True.
      else
         Is_First_Call = .False.
      end if
   end function Is_First_Call
