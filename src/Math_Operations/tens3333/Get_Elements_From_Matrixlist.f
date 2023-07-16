   ! --------------------------------------------------------------- !
   pure function Get_Elements_From_Matrixlist(matlist, nel, idx,jdx)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 3, 3), intent(in) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel) :: Get_Elements_From_Matrixlist
      ! ------------------------------------------------------------ !
      Get_Elements_From_Matrixlist = matlist(:, idx, jdx)
   end function Get_Elements_From_Matrixlist
