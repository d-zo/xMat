   ! --------------------------------------------------------------- !
   pure subroutine LU_Solve(lu_mat, nel, trans_mat, b_vec)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: lu_mat            ! Use LU decomposition of matrix `\mathbf{A}` saved in lu_mat and trans_mat
      real(dp), dimension(nel, nel), intent(in) :: trans_mat         ! representing `\mathbf{P}\mathbf{A} = \mathbf{L}\mathbf{U}` to calculate `\mathbf{A}x = b` with `\mathbf{L}y = \mathbf{P}b` and `\mathbf{U}x = y`.
      real(dp), dimension(nel), intent(inout) :: b_vec
      ! ------------------------------------------------------------ !
      integer :: idx

      ! All results are directly be written to b_mat, but association is used for code clarity
      associate(pb_vec => b_vec, y_vec => b_vec, x_vec => b_vec)
         pb_vec = matmul(trans_mat, pb_vec)                          ! `b^\prime = \mathbf{P}b`

         do idx = 2, nel                                             ! Solve for `y` in `\mathbf{L}y = \mathbf{P}b` by forward substitution
            y_vec(idx) = pb_vec(idx) &                               ! `y_i = \left(b_i^\prime \sum_{j=1}^{i-1} L_{ij}y_j\right) \frac{1}{L_{ii}}` with `L_{ii} = 1` and therefore `y_1 = b_1^\prime`
                       - dot_product(lu_mat(idx, 1:idx-1), y_vec(1:idx-1))
         end do

         x_vec(nel) = y_vec(nel)/lu_mat(nel, nel)                    ! Solve for `x` in `\mathbf{U}x = y` by back substitution
         do idx = nel-1, 1, -1                                       ! `x_i = \left(y_i - \sum_{j=i+1}^{n} U_{ij} x_j\right) \frac{1}{U_{ii}}`
            x_vec(idx) = (y_vec(idx) - dot_product(lu_mat(idx, idx+1:nel), x_vec(idx+1:nel)))/lu_mat(idx, idx)
         end do
         ! By association, all values determined for x_vec are already stored in b_vec
      end associate
   end subroutine LU_Solve
