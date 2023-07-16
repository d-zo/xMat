   ! --------------------------------------------------------------- !
   pure function Get_Elements_From_Matrixlist(matlist, nel, idx, jdx)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 6), intent(in) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel) :: Get_Elements_From_Matrixlist
      ! ------------------------------------------------------------ !
      Get_Elements_From_Matrixlist = matlist(:, ref_elements(idx, jdx))
   end function Get_Elements_From_Matrixlist
