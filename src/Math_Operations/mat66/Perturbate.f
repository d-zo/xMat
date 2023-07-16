   ! --------------------------------------------------------------- !
   pure subroutine Perturbate(mat, idx, jdx, theta)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(inout) :: mat
      integer, intent(in) :: idx, jdx
      real(dp), intent(in) :: theta
      ! ------------------------------------------------------------ !
      integer :: comb_idx

      comb_idx = ref_elements(idx, jdx)
      if (idx == jdx) then
         mat(comb_idx) = mat(comb_idx) + theta
      else
         mat(comb_idx) = mat(comb_idx) + 0.5_dp*theta
      end if
   end subroutine Perturbate
