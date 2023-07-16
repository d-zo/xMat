   ! --------------------------------------------------------------- !
   pure function Dyadic_Product22(mat33a, mat33b) result(tens_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33a, mat33b
      real(dp), dimension(3, 3, 3, 3) :: tens_out
      ! ------------------------------------------------------------ !
      integer :: kdx, ldx

      do ldx = 1, 3
         do kdx = 1, 3
            tens_out(:, :, kdx, ldx) = mat33a*mat33b(kdx, ldx)       ! Dyadic product `T_{ijkl} = M_{A,ij}\cdot M_{B,kl}`
         end do
      end do
   end function Dyadic_Product22
