   ! --------------------------------------------------------------- !
   pure subroutine Set_Element_In_Tensor(tens, idx, jdx, kdx, ldx, val)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(inout) :: tens
      integer, intent(in) :: idx, jdx, kdx, ldx
      real(dp), intent(in) :: val
      ! ------------------------------------------------------------ !
      tens(idx, jdx, kdx, ldx) = val
   end subroutine Set_Element_In_Tensor
