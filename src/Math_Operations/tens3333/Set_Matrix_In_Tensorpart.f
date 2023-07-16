   ! --------------------------------------------------------------- !
   pure subroutine Set_Matrix_In_Tensorpart(tens, idx, jdx, mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(inout) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(3, 3), intent(in) :: mat
      ! ------------------------------------------------------------ !
      tens(:, :, idx, jdx) = mat
   end subroutine Set_Matrix_In_Tensorpart
