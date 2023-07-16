   ! --------------------------------------------------------------- !
   pure subroutine Set_Elements_In_Matrixlist(matlist, nel, idx, jdx, elemlist)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 3, 3), intent(inout) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel), intent(in) :: elemlist
      ! ------------------------------------------------------------ !
      matlist(:, idx, jdx) = elemlist
   end subroutine Set_Elements_In_Matrixlist
