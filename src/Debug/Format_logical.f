   ! --------------------------------------------------------------- !
   pure function Format_logical(name, bool)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      logical, intent(in) :: bool
      character(len=len(name)+6) :: Format_logical
      ! ------------------------------------------------------------ !
      if (bool) then
         write(Format_logical, '(a, a)') name, ':  Yes'
      else
         write(Format_logical, '(a, a)') name, ':  No '
      end if
   end function Format_logical
