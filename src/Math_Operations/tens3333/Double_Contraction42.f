   ! --------------------------------------------------------------- !
   pure function Double_Contraction42(tens3333, mat33) result(mat_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tens3333
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp), dimension(3, 3) :: mat_out
      ! ------------------------------------------------------------ !
      integer :: idx, jdx

      mat_out = 0.0_dp                                               ! Double contraction `\mathcal{T} : \mathbf{M} = \sum_{kl} T_{ijkl}\cdot M_{kl}`
      do idx = 1, 3
         do jdx = 1, 3
            mat_out = mat_out + reshape(tens3333(:, :, jdx, idx), [3, 3])*mat33(jdx, idx)
         end do
      end do
   end function Double_Contraction42
