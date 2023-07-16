   ! --------------------------------------------------------------- !
   pure function Uppercase(text) result(upper)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: text
      character(len=len(text)) :: upper
      ! ------------------------------------------------------------ !
      character(len=26), parameter :: low_alph = 'abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: up_alph  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer :: idx, pos

      upper = text
      do idx = 1, len(text)
         pos = index(low_alph, text(idx:idx))
         if (pos /= 0) then
            upper(idx:idx) = up_alph(pos:pos)
         end if
      end do
   end function Uppercase
