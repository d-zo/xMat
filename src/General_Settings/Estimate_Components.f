   ! --------------------------------------------------------------- !
   pure subroutine Estimate_Components(nel, num_dimensions, num_shear)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel                                     ! The number of elements nel has to be between 3 and 6
      integer, intent(out) :: num_dimensions, num_shear
      ! ------------------------------------------------------------ !
      num_dimensions = 0
      num_shear = 0

      if (nel == 3) then                                             ! Either two direct and one shear component
         num_dimensions = 2
         num_shear = 1
      else if ((nel > 3) .and. (nel <= 6)) then                      ! Or three direct and up to three shear components for 3 < nel <= 6
         num_dimensions = 3
         num_shear = nel - num_dimensions
      end if
   end subroutine Estimate_Components
