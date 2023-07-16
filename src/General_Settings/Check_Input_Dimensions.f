   ! --------------------------------------------------------------- !
   function Check_Input_Dimensions(num_dimensions, num_shear)
   ! --------------------------------------------------------------- ! Only support 2D and 3D symmetric input
      integer, intent(in) :: num_dimensions, num_shear
      logical :: Check_Input_Dimensions
      ! ------------------------------------------------------------ !
      Check_Input_Dimensions = .False.

      if (num_dimensions == 2) then
         if (num_shear == 1) then
            Check_Input_Dimensions = .True.
         end if
      else if (num_dimensions == 3) then
         if ((num_shear >= 1) .and. (num_shear <= 3)) then
            Check_Input_Dimensions = .True.
         end if
      end if

      if (Check_Input_Dimensions) then
         global_num_direct_components = num_dimensions               ! Changing global variables makes this function impure
         global_num_shear_components = num_shear
      end if
   end function Check_Input_Dimensions
