   ! --------------------------------------------------------------- !
   pure function Get_Matrix_From_Tensorpart(tens, idx, jdx)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(9) :: Get_Matrix_From_Tensorpart
      ! ------------------------------------------------------------ !
      Get_Matrix_From_Tensorpart = tens(:, ref_elements(idx, jdx))
   end function Get_Matrix_From_Tensorpart
