   ! --------------------------------------------------------------- !
   pure function Format_int_dim0(name, number)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      integer, intent(in) :: number
      character(len=len(name)+2 + 6) :: Format_int_dim0
      ! ------------------------------------------------------------ !
      write(Format_int_dim0, '(a, a, i6)') name, '  ', number
   end function Format_int_dim0
