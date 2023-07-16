   ! --------------------------------------------------------------- !
   pure function Is_Diagonal(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      logical :: Is_Diagonal
      ! ------------------------------------------------------------ !
      integer :: idx, jdx

      Is_Diagonal = .True.
      do jdx = 1, 3
         do idx = 1, 3
            if ((idx /= jdx) .and. (abs(mat33(idx, jdx)) > setting_epsilon)) then
               Is_Diagonal = .False.
            end if
         end do
      end do
   end function Is_Diagonal
