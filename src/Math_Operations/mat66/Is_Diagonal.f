   ! --------------------------------------------------------------- !
   pure function Is_Diagonal(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      logical :: Is_Diagonal
      ! ------------------------------------------------------------ !
      if (any(abs(vec6(4:6)) > setting_epsilon)) then
         Is_Diagonal = .False.
      else
         Is_Diagonal = .True.
      end if
   end function Is_Diagonal
