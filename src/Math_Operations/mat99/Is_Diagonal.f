   ! --------------------------------------------------------------- !
   pure function Is_Diagonal(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      logical :: Is_Diagonal
      ! ------------------------------------------------------------ !
      if (any(abs(vec9(4:9)) > setting_epsilon)) then
         Is_Diagonal = .False.
      else
         Is_Diagonal = .True.
      end if
   end function Is_Diagonal
