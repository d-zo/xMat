   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Internal(matrix, inv_matrix, nel, success)
   ! --------------------------------------------------------------- ! Calculates the inverse of a nel by nel matrix if it is not singular
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: matrix
      real(dp), dimension(nel, nel), intent(out) :: inv_matrix
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: lu_mat, p_mat
      integer :: idx, jdx

      lu_mat = matrix
      inv_matrix = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])

      ! Calculates LU-decomposition of `\mathbf{A}` and solves `Ax_j = \delta_j` for each column vector `\delta_j` of the identity matrix to get `\left[x_1| \dots |x_n\right] = \mathbf{A}^{-1}`
      call LU_Decomposition(matrix=lu_mat, nel=size(matrix, 1), trans_mat=p_mat, success=success)
      if (success) then
         do idx = 1, nel
            call LU_Solve(lu_mat=lu_mat, nel=nel, trans_mat=p_mat, b_vec=inv_matrix(:, idx))
         end do
      end if
   end subroutine Inverse_Internal
