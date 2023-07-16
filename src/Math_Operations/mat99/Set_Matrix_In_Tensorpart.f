   ! --------------------------------------------------------------- !
   pure subroutine Set_Matrix_In_Tensorpart(tens, idx, jdx, mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(inout) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(9), intent(in) :: mat
      ! ------------------------------------------------------------ !
      tens(:, ref_elements(idx, jdx)) = mat
   end subroutine Set_Matrix_In_Tensorpart
