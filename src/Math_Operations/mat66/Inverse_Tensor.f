   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Tensor(tensor, inv_tensor, successful_inversion)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: tensor
      real(dp), dimension(6, 6), intent(out) :: inv_tensor
      logical, intent(out) :: successful_inversion
      ! ------------------------------------------------------------ !
      real(dp), dimension(global_num_direct_components+global_num_shear_components, &
         global_num_direct_components+global_num_shear_components) :: comp_matrix, inv_comp_matrix
      logical :: is_singular

      ! To successfully find an inverse, the matrix has to have full rank. Since all calculations are done with vec6 or mat66,
      ! for inversion all rows and columns added for simplicity have to be temporarily removed.
      ! If the inversion fails, successful_inversion is set to false and the calling function has to deal with it.
      associate(ndi => global_num_direct_components, nshr => global_num_shear_components)
         comp_matrix(1:ndi, 1:ndi) = tensor(1:ndi, 1:ndi)
         comp_matrix((ndi + 1):(ndi + nshr),1:ndi) = tensor(4:(3 + nshr), 1:ndi)
         comp_matrix(1:ndi, (ndi + 1):(ndi + nshr)) = tensor(1:ndi, 4:(3 + nshr))
         comp_matrix((ndi + 1):(ndi + nshr), (ndi + 1):(ndi + nshr)) = tensor(4:(3 + nshr), 4:(3 + nshr))

         call Inverse_Internal(matrix=comp_matrix, inv_matrix=inv_comp_matrix, nel=ndi+nshr, is_singular=is_singular)
         successful_inversion = .not. is_singular                    ! Only regular matrices have a successful inversion

         inv_tensor = 0.0_dp
         inv_tensor(1:ndi, 1:ndi) = inv_comp_matrix(1:ndi, 1:ndi)
         inv_tensor(4:(3 + nshr),1:ndi) = inv_comp_matrix((ndi + 1):(ndi + nshr), 1:ndi)
         inv_tensor(1:ndi, 4:(3 + nshr)) = inv_comp_matrix(1:ndi, (ndi + 1):(ndi + nshr))
         inv_tensor(4:(3 + nshr), 4:(3 + nshr)) = inv_comp_matrix((ndi + 1):(ndi + nshr), (ndi + 1):(ndi + nshr))
      end associate
   end subroutine Inverse_Tensor
