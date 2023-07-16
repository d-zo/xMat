   ! --------------------------------------------------------------- !
   pure subroutine Ref_Index(ref_idx, idx, jdx)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: ref_idx
      integer, intent(out) :: idx, jdx
      ! ------------------------------------------------------------ !
      idx = ref_indices(1, ref_idx)
      jdx = ref_indices(2, ref_idx)
   end subroutine Ref_Index
