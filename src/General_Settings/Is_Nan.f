   ! --------------------------------------------------------------- !
   elemental function Is_Nan(number)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: number
      logical :: Is_Nan
      ! ------------------------------------------------------------ !
      if (number /= number) then                                     ! Only true for NaN
         Is_Nan = .True.
      else
         Is_Nan = .False.
      end if
   end function Is_Nan
