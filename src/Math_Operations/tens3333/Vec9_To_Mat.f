   ! --------------------------------------------------------------- !
   pure function Vec9_To_Mat(vec9) result(mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(3, 3) :: mat
      ! ------------------------------------------------------------ !
      integer :: ref_idx, idx, jdx

      do ref_idx = 1, 9
         call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
         mat(idx, jdx) = vec9(ref_idx)
      end do
   end function Vec9_To_Mat
