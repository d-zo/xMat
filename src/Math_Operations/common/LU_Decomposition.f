   ! --------------------------------------------------------------- !
   pure subroutine LU_Decomposition(matrix, nel, trans_mat, success)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! Determine the LU-decomposition of matrix `\mathbf{A}`, so that `\mathbf{P}\mathbf{A} = \mathbf{L}\mathbf{U}`
      real(dp), dimension(nel, nel), intent(out) :: trans_mat        ! Save `\mathbf{L}` and `\mathbf{U}` in matrix and the permutation `\mathbf{P}` in trans_mat
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), parameter, dimension(2, 2) :: permutation = reshape([ &
         0.0_dp, 1.0_dp, &
         1.0_dp, 0.0_dp], [2, 2])
      integer :: idx, jdx, idx_max

      success = .True.
      trans_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])

      do idx = 1, nel-1
         ! Find the largest element for (partial) pivoting and adjust transformation matrix
         idx_max = idx + maxloc(abs(matrix(idx:nel, idx)), dim=1) - 1
         if (idx_max /= idx) then
            trans_mat([idx, idx_max], :) = matmul(permutation, trans_mat([idx, idx_max], :))
            matrix([idx, idx_max], :) = matrix([idx_max, idx], :)
         end if

         if (abs(matrix(idx, idx)) < setting_epsilon_extra) then
            ! If value on diagonal (maximum value of matrix(idx:nel, idx)) is zero, the matrix is singular
            success = .False.
            cycle
         end if

         matrix(idx+1:nel, idx) = matrix(idx+1:nel, idx)/matrix(idx, idx)
         matrix(idx+1:nel, idx+1:nel) = matrix(idx+1:nel, idx+1:nel) &
                                      - matmul(reshape(matrix(idx+1:nel, idx), [nel-idx, 1]), &
                                               reshape(matrix(idx, idx+1:nel), [1, nel-idx]))
      end do
   end subroutine LU_Decomposition
