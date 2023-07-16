   ! --------------------------------------------------------------- !
   elemental function Number_To_String(number)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan
      !
      real(dp), intent(in) :: number
      character(len=onenumberlength+2) :: Number_To_String
      ! ------------------------------------------------------------ !
      character(len=len(debug_format)+2) :: formatstring

      Number_To_String(1:onenumberlength+2) = repeat(' ', onenumberlength+2)
      Number_To_String(onenumberlength+1:onenumberlength+1) = ','

      if (Is_Nan(number)) then
         write(Number_To_String(onenumberlength-2:onenumberlength), '(a)') 'NaN'
      else                                                           ! Use fixed format for numbers "close" to zero
         if ((abs(number) > 1000.0_dp) .or. (abs(number) < 0.001_dp)) then
            if (abs(number) < setting_epsilon) then
               write(Number_To_String(onenumberlength-2:onenumberlength), '(a)') '0.0'
            else
               write(formatstring, '(a, a, a)') '(', alternative_format, ')'
               write(Number_To_String(1:onenumberlength), fmt=formatstring) number
            end if
         else
            write(formatstring, '(a, a, a)') '(', debug_format, ')'
            write(Number_To_String(1:onenumberlength), fmt=formatstring) number
         end if
      end if
   end function Number_To_String
