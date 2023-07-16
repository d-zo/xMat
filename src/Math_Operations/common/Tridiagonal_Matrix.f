   ! --------------------------------------------------------------- !
   pure subroutine Tridiagonal_Matrix(matrix, nel, trans_mat)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! Computes a tridiagonal matrix and an assembly of transformation
      real(dp), dimension(nel, nel), intent(out) :: trans_mat        ! matrices as described in section 11.2 of Press et al (1997)
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: ident_mat, householder_mat    ! The explicit creation of the householder matrices (as done here)
      real(dp), dimension(nel) :: u_vec                              ! is only justifiable for small matrices (nel << 10)
      real(dp) :: H_mat, sigma, scaling_factor
      integer :: idx, jdx, ldx

      ident_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])
      trans_mat = ident_mat                                          ! `\mathbf{Q}_0 = \mathbf{I}`

      do idx = nel, 3, -1
         ldx = idx - 1
         scaling_factor = sum(abs(matrix(idx, 1:ldx)))

         if (scaling_factor > 0.0_dp) then
            sigma = dot_product(matrix(idx, 1:ldx)/scaling_factor, matrix(idx, 1:ldx)/scaling_factor)
            sigma = sign(sqrt(sigma), matrix(idx, ldx))

            u_vec = 0.0_dp
            u_vec(1:ldx) = matrix(idx, 1:ldx)/scaling_factor
            u_vec(ldx) = u_vec(ldx) + sigma
            H_mat = 0.5_dp*dot_product(u_vec(1:ldx), u_vec(1:ldx))   ! (11.2.4) of Press et al (1997): `H = \frac{1}{2}|u|^2`

            householder_mat = ident_mat                              ! (11.2.3) of Press et al (1997): `\mathbf{P}_{n+1} = \mathbf{I} - \frac{u u^T}{H}`
            householder_mat(1:ldx, 1:ldx) = householder_mat(1:ldx, 1:ldx) &
                                          - matmul(reshape(u_vec(1:ldx)/H_mat, [ldx, 1]), &
                                                   reshape(u_vec(1:ldx), [1, ldx]))
            matrix = matmul(householder_mat, matmul(matrix, transpose(householder_mat)))
            trans_mat = matmul(trans_mat, householder_mat)           ! (11.2.17) of Press et al (1997): `\mathbf{Q}_{n+1} = \mathbf{Q}_n \mathbf{P}_{n+1}`
         end if
      end do

      ! Set all elements above/beneath subdiagonal to (exactly) zero
      do idx = 1, nel
         do jdx = 1, nel
            if (abs(jdx-idx) > 1) then
               matrix(jdx, idx) = 0.0_dp
            end if
         end do
      end do
   end subroutine Tridiagonal_Matrix
