   ! --------------------------------------------------------------- !
   pure subroutine Set_Element_In_Tensor(tens, idx, jdx, kdx, ldx, val)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(inout) :: tens
      integer, intent(in) :: idx, jdx, kdx, ldx
      real(dp), intent(in) :: val
      ! ------------------------------------------------------------ !
      tens(ref_elements(idx, jdx), ref_elements(kdx, ldx)) = val
   end subroutine Set_Element_In_Tensor
