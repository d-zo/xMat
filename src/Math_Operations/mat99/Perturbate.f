   ! --------------------------------------------------------------- !
   pure subroutine Perturbate(mat, idx, jdx, theta)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(inout) :: mat
      integer, intent(in) :: idx, jdx
      real(dp), intent(in) :: theta
      ! ------------------------------------------------------------ !
      mat(ref_elements(idx, jdx)) = mat(ref_elements(idx, jdx)) + theta
   end subroutine Perturbate
