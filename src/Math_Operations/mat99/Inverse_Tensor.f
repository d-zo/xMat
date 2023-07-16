   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Tensor(tensor, inv_tensor, successful_inversion)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9, 9), intent(in) :: tensor
      real(dp), dimension(9, 9), intent(out) :: inv_tensor
      logical, intent(out) :: successful_inversion
      ! ------------------------------------------------------------ !
      real(dp), dimension(global_num_direct_components+global_num_shear_components, &
         global_num_direct_components+global_num_shear_components) :: comp_matrix, inv_comp_matrix
      integer :: idx, jdx, ref_jdx
      logical :: is_singular

      ! To successfully find an inverse, the matrix has to have full rank. If the amount of direct and shear components
      ! is less than six, additional rows and columns used before for simplicity have to be temporarily removed.
      ! If the inversion fails, successful_inversion is set to false and the calling function has to deal with it.
      associate(ndi => global_num_direct_components, nshr => global_num_shear_components)
         comp_matrix(1:ndi, 1:ndi) = tensor(1:ndi, 1:ndi)
         do jdx = 1, nshr
            ref_jdx = 2+2*jdx
            do idx = 1, ndi
               comp_matrix(idx, 3+jdx) = 0.5_dp*(tensor(idx, ref_jdx) + tensor(idx, ref_jdx+1))
               comp_matrix(3+jdx, idx) = 0.5_dp*(tensor(ref_jdx, idx) + tensor(ref_jdx+1, idx))
            end do
            comp_matrix(3+jdx, 3+jdx) = 0.5_dp*(tensor(ref_jdx, ref_jdx) + tensor(ref_jdx, ref_jdx+1) &
                                      + tensor(ref_jdx+1, ref_jdx) + tensor(ref_jdx+1, ref_jdx+1))
         end do

         call Inverse_Internal(matrix=comp_matrix, inv_matrix=inv_comp_matrix, nel=ndi+nshr, is_singular=is_singular)
         successful_inversion = .not. is_singular                    ! Only regular matrices have a successful inversion

         inv_tensor = 0.0_dp
         inv_tensor(1:ndi, 1:ndi) = inv_comp_matrix(1:ndi, 1:ndi)
         do jdx = 1, nshr
            ref_jdx = 2+2*jdx
            do idx = 1, ndi
               inv_tensor(idx, ref_jdx) = inv_comp_matrix(idx, 3+jdx)
               inv_tensor(idx, ref_jdx+1) = inv_comp_matrix(idx, 3+jdx)
               inv_tensor(ref_jdx, idx) = inv_comp_matrix(3+jdx, idx)
               inv_tensor(ref_jdx+1, idx) = inv_comp_matrix(3+jdx, idx)
            end do
            inv_tensor(ref_jdx, ref_jdx) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_tensor(ref_jdx+1, ref_jdx) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_tensor(ref_jdx, ref_jdx+1) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_tensor(ref_jdx+1, ref_jdx+1) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
         end do
      end associate
   end subroutine Inverse_Tensor
