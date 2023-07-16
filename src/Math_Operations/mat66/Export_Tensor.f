   ! --------------------------------------------------------------- !
   pure function Export_Tensor(tens, num_dimensions, num_shear) result(tens_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of tens_out may need to be adjusted hereafter
      real(dp), dimension(6, 6), intent(in) :: tens
      real(dp), dimension(num_dimensions+num_shear, num_dimensions+num_shear) :: tens_out
      ! ------------------------------------------------------------ !
      integer :: nel

      nel = num_dimensions + num_shear
      tens_out = 0.0_dp
      tens_out(1:num_dimensions, 1:num_dimensions) = tens(1:num_dimensions, 1:num_dimensions)
      tens_out((num_dimensions + 1):nel, 1:num_dimensions) = tens(4:(3 + num_shear), 1:num_dimensions)
      tens_out(1:num_dimensions, (num_dimensions + 1):nel) = tens(1:num_dimensions, 4:(3 + num_shear))
      tens_out((num_dimensions + 1):nel, (num_dimensions + 1):nel) = tens(4:(3 + num_shear), 4:(3 + num_shear))
   end function Export_Tensor
