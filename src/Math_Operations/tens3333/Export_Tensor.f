   ! --------------------------------------------------------------- !
   pure function Export_Tensor(tens, num_dimensions, num_shear) result(tens_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of tens_out may need to be adjusted hereafter
      real(dp), dimension(3, 3, 3, 3), intent(in) :: tens
      real(dp), dimension(num_dimensions+num_shear, num_dimensions+num_shear) :: tens_out
      ! ------------------------------------------------------------ !
      integer, dimension(4), parameter :: indexmap = [1, 2, 3, 1]
      integer :: idx, jdx, nel

      associate (ndi => num_dimensions, nshr => num_shear)
         nel = ndi + nshr
         tens_out = 0.0_dp
         do jdx = 1, ndi
            do idx = 1, ndi
               tens_out(idx, jdx) = tens(idx, idx, jdx, jdx)
            end do

            do idx = 1, nshr
               tens_out(ndi+idx, jdx) = 0.5_dp*(tens(indexmap(idx), indexmap(idx+1), jdx, jdx) &
                                                + tens(indexmap(idx+1), indexmap(idx), jdx, jdx))
               tens_out(jdx, ndi+idx) = 0.5_dp*(tens(jdx, jdx, indexmap(idx), indexmap(idx+1)) &
                                                + tens(jdx, jdx, indexmap(idx+1), indexmap(idx)))
            end do
         end do

         do jdx = 1, nshr
            do idx = 1, nshr
               tens_out(ndi+idx, ndi+jdx) = 0.25_dp*( &
                  tens(indexmap(idx), indexmap(idx+1), indexmap(jdx), indexmap(jdx+1)) &
                  + tens(indexmap(idx), indexmap(idx+1), indexmap(jdx+1), indexmap(jdx)) &
                  + tens(indexmap(idx+1), indexmap(idx), indexmap(jdx), indexmap(jdx+1)) &
                  + tens(indexmap(idx+1), indexmap(idx), indexmap(jdx+1), indexmap(jdx)))
            end do
         end do
      end associate
   end function Export_Tensor
