   ! --------------------------------------------------------------- !
   function Matrix_Exponential(mat) result(mat_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat
      real(dp), dimension(3, 3) :: mat_out
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: temp_mat, eigenval_diag, eigenvec
      integer :: idx

      ! `e^{a\mathbf{M}}` is meant to apply the `e^ax`-function to each entry `x` of the diagonal matrix `\mathbf{M}`
      ! If `\mathbf{M}` is not diagonal, decompose its eigenvalues/eigenvectors first, so that `e^{a\mathbf{M}} = \mathbf{V}\cdot\mathrm{diag}\left(e^{a\lambda_1},\dots,e^{a\lambda_n}\right)\mathbf{V}^{-1}`
      if (Is_Diagonal(mat)) then
         mat_out = 0.0_dp
         do idx = 1, 3
            mat_out(idx, idx) = exp(mat(idx, idx))
         end do
      else
         call Eigendecomposition(mat=mat, nel=3, eigenvalues=eigenval_diag, eigenvectors=eigenvec)

         temp_mat = 0.0_dp
         do idx = 1, 3
            temp_mat(idx, idx) = exp(eigenval_diag(idx, idx))
         end do
         mat_out = matmul(eigenvec, matmul(temp_mat, transpose(eigenvec)))
      end if
   end function Matrix_Exponential
