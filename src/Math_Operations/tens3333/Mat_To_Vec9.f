   ! --------------------------------------------------------------- !
   pure function Mat_To_Vec9(mat) result(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat
      real(dp), dimension(9) :: vec9
      ! ------------------------------------------------------------ !
      integer :: ref_idx, idx, jdx

      do ref_idx = 1, 9
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         vec9(ref_idx) = mat(idx, jdx)
      end do
   end function Mat_To_Vec9
