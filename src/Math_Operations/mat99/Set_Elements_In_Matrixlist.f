   ! --------------------------------------------------------------- !
   pure subroutine Set_Elements_In_Matrixlist(matlist, nel, idx, jdx, elemlist)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 9), intent(inout) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel), intent(in) :: elemlist
      ! ------------------------------------------------------------ !
      matlist(:, ref_elements(idx, jdx)) = elemlist
   end subroutine Set_Elements_In_Matrixlist
