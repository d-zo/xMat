   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Tensor(tensor, inv_tensor, success)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tensor
      real(dp), dimension(3, 3, 3, 3), intent(out) :: inv_tensor
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), dimension(9, 9) :: proj_matrix
      real(dp), dimension(9, 9) :: inv_proj_matrix
      real(dp), dimension(global_num_direct_components+global_num_shear_components, &
         global_num_direct_components+global_num_shear_components) :: comp_matrix, inv_comp_matrix
      integer :: idx, jdx, ref_jdx

      ! To successfully find an inverse, the matrix has to have full rank. If the amount of direct and shear components
      ! is less than six, additional rows and columns used before for simplicity have to be temporarily removed.
      ! If the inversion fails, success is set to .False. and the calling function has to deal with it.
      associate(ndi => global_num_direct_components, nshr => global_num_shear_components)
         proj_matrix = Tens_To_Mat99(tensor)
         comp_matrix(1:ndi, 1:ndi) = proj_matrix(1:ndi, 1:ndi)
         do jdx = 1, nshr
            ref_jdx = 2+2*jdx
            do idx = 1, ndi
               comp_matrix(idx, 3+jdx) = 0.5_dp*(proj_matrix(idx, ref_jdx) + proj_matrix(idx, ref_jdx+1))
               comp_matrix(3+jdx, idx) = 0.5_dp*(proj_matrix(ref_jdx, idx) + proj_matrix(ref_jdx+1, idx))
            end do
            comp_matrix(3+jdx, 3+jdx) = 0.5_dp*(proj_matrix(ref_jdx, ref_jdx) + proj_matrix(ref_jdx, ref_jdx+1) &
                                      + proj_matrix(ref_jdx+1, ref_jdx) + proj_matrix(ref_jdx+1, ref_jdx+1))
         end do

         call Inverse_Internal(matrix=comp_matrix, inv_matrix=inv_comp_matrix, nel=ndi+nshr, success=success)

         inv_proj_matrix = 0.0_dp
         inv_proj_matrix(1:ndi, 1:ndi) = inv_comp_matrix(1:ndi, 1:ndi)
         do jdx = 1, nshr
            ref_jdx = 2+2*jdx
            do idx = 1, ndi
               inv_proj_matrix(idx, ref_jdx) = inv_comp_matrix(idx, 3+jdx)
               inv_proj_matrix(idx, ref_jdx+1) = inv_comp_matrix(idx, 3+jdx)
               inv_proj_matrix(ref_jdx, idx) = inv_comp_matrix(3+jdx, idx)
               inv_proj_matrix(ref_jdx+1, idx) = inv_comp_matrix(3+jdx, idx)
            end do
            inv_proj_matrix(ref_jdx, ref_jdx) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_proj_matrix(ref_jdx+1, ref_jdx) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_proj_matrix(ref_jdx, ref_jdx+1) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
            inv_proj_matrix(ref_jdx+1, ref_jdx+1) = 0.5_dp*inv_comp_matrix(3+jdx, 3+jdx)
         end do
         inv_tensor = Mat99_To_Tens(inv_proj_matrix)
      end associate
   end subroutine Inverse_Tensor
