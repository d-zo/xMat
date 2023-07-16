   ! --------------------------------------------------------------- !
   pure function C_To_Fortran_String(length, c_string) result(f_string)
   ! --------------------------------------------------------------- !
      use iso_c_binding, only: c_char, c_null_char

      integer, intent(in) :: length
      character(kind=c_char, len=1), dimension(length), intent(in) :: c_string
      character(len=length) :: f_string
      ! ------------------------------------------------------------ !
      integer :: idx

      f_string = " "
      do idx = 1, length
         if (c_string(idx) == c_null_char) then
            exit
         else
            f_string(idx:idx) = c_string(idx)
         end if
      end do
   end function C_To_Fortran_String
