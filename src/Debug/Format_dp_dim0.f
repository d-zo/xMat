   ! --------------------------------------------------------------- !
   pure function Format_dp_dim0(name, number)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: number
      character(len=len(name)+2 + onenumberlength+2) :: Format_dp_dim0
      ! ------------------------------------------------------------ !
      write(Format_dp_dim0, '(a, a, a)') name, '  ', Number_To_String(number)
      Format_dp_dim0(len(Format_dp_dim0)-1:len(Format_dp_dim0)-1) = ' '
   end function Format_dp_dim0
