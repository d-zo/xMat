   ! --------------------------------------------------------------- !
   pure subroutine Perturbate(mat, idx, jdx, theta)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(inout) :: mat
      integer, intent(in) :: idx, jdx
      real(dp), intent(in) :: theta
      ! ------------------------------------------------------------ !
      mat(idx, jdx) = mat(idx, jdx) + theta
   end subroutine Perturbate
