   ! --------------------------------------------------------------- !
   pure function Get_Matrix_From_Tensorpart(tens, idx, jdx)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(3, 3) :: Get_Matrix_From_Tensorpart
      ! ------------------------------------------------------------ !
      Get_Matrix_From_Tensorpart = tens(:, :, idx, jdx)
   end function Get_Matrix_From_Tensorpart
