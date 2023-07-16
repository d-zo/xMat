   ! --------------------------------------------------------------- !
   pure function Format_dp_dim1(name, vector)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:), intent(in) :: vector
      character(len=len(name)+2 + (onenumberlength+2)*size(vector)) :: Format_dp_dim1
      ! ------------------------------------------------------------ !
      character(len=11) :: formatstring

      write(formatstring, '(a, i2, a)') '(a, a, ', size(vector), 'a)'
      write(Format_dp_dim1, formatstring) name, '  ', Number_To_String(vector)
      Format_dp_dim1(len(Format_dp_dim1)-1:len(Format_dp_dim1)-1) = ' '
      Format_dp_dim1(len(Format_dp_dim1):len(Format_dp_dim1)) = char(10)
   end function Format_dp_dim1
