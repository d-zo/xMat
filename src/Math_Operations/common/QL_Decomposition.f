   ! --------------------------------------------------------------- !
   pure subroutine QL_Decomposition(matrix, nel, eigenvec_mat, success)
   ! --------------------------------------------------------------- ! Use this function for symmetric tridiagonal matrices only.
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! QL decomposition with implicit shifts of matrix as described in
      real(dp), dimension(nel, nel), intent(out) :: eigenvec_mat     ! section 11.3 of Press et al (1997) returning the eigenvalues in
      logical, intent(out) :: success                                ! matrix (not sorted), the eigenvectors and a success flag
      ! ------------------------------------------------------------ !
      integer, parameter :: max_iterations = 30
      integer :: idx, jdx, kdx, idx_start, idx_end, iter
      real(dp) :: diff, mu, refval, zeroval, temp_sum
      real(dp), dimension(2, 2) :: givens
      real(dp), parameter, dimension(2, 2) :: identity_2d = reshape([ &
         1.0_dp, 0.0_dp, &
         0.0_dp, 1.0_dp], [2, 2])
      logical :: extract_block

      eigenvec_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])
      success = .False.

      do idx_start = 1, nel
         iter = 0
         do
            iter = iter+1
            if (iter == max_iterations) then
               ! Too many iterations: Cancel the decomposition procedure
               return
            end if

            extract_block = .False.
            ! Look for a single small subdiagonal element to split the matrix in two blocks.
            ! Focus on the top left block starting from (idx_start, idx_start)
            do idx_end = idx_start, nel-1
               temp_sum = abs(matrix(idx_end, idx_end)) + abs(matrix(idx_end+1, idx_end+1))
               ! If adding the subdiagonal element to temp_sum is smaller than the next representable
               ! value in the used precision, then it is negligible and the matrix can be split here
               if (nearest(temp_sum, 1.0_dp) >= temp_sum + abs(matrix(idx_end+1, idx_end))) then
                  extract_block = .True.
                  exit
               end if
            end do

            if (.not. extract_block) then
               idx_end = nel
            end if

            if (idx_end == idx_start) then
               exit
            end if

            ! The first rotation is special while following (idx_end-1-idx_start) rotations "chase"
            ! the only non-zero off-subdiagonal element out of the matrix
            diff = (matrix(idx_start+1, idx_start+1) - matrix(idx_start, idx_start))/2.0_dp
            mu = matrix(idx_start, idx_start) - matrix(idx_start+1, idx_start)**2 &
               / (diff + sign(Sqrt_Of_Sum_Of_Squares(diff, matrix(idx_start+1, idx_start)), diff))
            refval = matrix(idx_end, idx_end) - mu
            zeroval = matrix(idx_end, idx_end-1)

            do kdx = idx_end-1, idx_start, -1
               if (kdx < idx_end-1) then
                  refval = matrix(kdx+2, kdx+1)
                  zeroval = matrix(kdx+2, kdx)
               end if

               givens = Givens_CS(refval=refval, zeroval=zeroval)
               if (kdx == idx_end-1) then
                  if (all(givens == identity_2d)) then
                     ! If the givens rotation is the same as the identity matrix in the used precision
                     ! then explicitly set the investigated value to zero and continue with the next value
                     matrix(idx_end, idx_end-1) = 0.0_dp
                     exit
                  end if
               end if

               ! Instead of creating the full matrices for a Givens rotation `\mathbf{G}\cdot\mathbf{A}\cdot\mathbf{G}^T`, only use the relevant
               ! slices (only row and column of kdx and kdx+1 change, but two commands seem necessary)
               matrix(:, [kdx, kdx+1]) = matmul(matrix(:, [kdx, kdx+1]), givens)
               matrix([kdx, kdx+1], :) = matmul(transpose(givens), matrix([kdx, kdx+1], :))

               eigenvec_mat(:, [kdx, kdx+1]) = matmul(eigenvec_mat(:, [kdx, kdx+1]), givens)
            end do
         end do
      end do

      ! All nondiagonal elements should be zero or close to it. Explicitly set them to zero and report success
      do idx = 1, nel
         refval = matrix(idx, idx)
         matrix(:, idx) = 0.0_dp
         matrix(idx, idx) = refval
      end do
      success = .True.
   end subroutine QL_Decomposition
