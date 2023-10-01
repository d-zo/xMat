   ! --------------------------------------------------------------- !
   subroutine Eigendecomposition(mat, nel, eigenvalues, eigenvectors)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: mat
      real(dp), dimension(nel, nel), intent(out) :: eigenvalues, eigenvectors
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: temp_mat, trans_mat
      logical :: success

      temp_mat = mat
      call Tridiagonal_Matrix(matrix=temp_mat, nel=nel, trans_mat=trans_mat)
      call QL_Decomposition(matrix=temp_mat, nel=nel, eigenvec_mat=eigenvectors, success=success)

      if (.not. success) then
         call Write_Error_And_Exit('Eigendecomposition: Failed QL decomposition of' // char(10) &
            // Formatval('mat', mat))
      end if

      eigenvectors = matmul(trans_mat, eigenvectors)
      eigenvalues = temp_mat
   end subroutine Eigendecomposition
