   ! --------------------------------------------------------------- !
   pure function Export_Tensor(tens, num_dimensions, num_shear) result(tens_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of tens_out may need to be adjusted hereafter
      real(dp), dimension(9, 9), intent(in) :: tens
      real(dp), dimension(num_dimensions+num_shear, num_dimensions+num_shear) :: tens_out
      ! ------------------------------------------------------------ !
      integer :: idx, jdx, nel

      associate (ndi => num_dimensions, nshr => num_shear)
         nel = ndi + nshr
         tens_out = 0.0_dp
         tens_out(1:ndi, 1:ndi) = tens(1:ndi, 1:ndi)

         do idx = 1, nshr
            tens_out(ndi+idx, 1:ndi) = 0.5_dp*(tens(2+2*idx, 1:ndi) + tens(3+2*idx, 1:ndi))
            tens_out(1:ndi, ndi+idx) = 0.5_dp*(tens(1:ndi, 2+2*idx) + tens(1:ndi, 3+2*idx))

            do jdx = 1, nshr
               tens_out(ndi+idx, ndi+jdx) = 0.25_dp*(tens(2+2*idx, 2+2*jdx) + tens(3+2*idx, 2+2*jdx) &
                                                     + tens(2+2*idx, 3+2*jdx) + tens(3+2*idx, 3+2*jdx))
            end do
         end do
      end associate
   end function Export_Tensor
