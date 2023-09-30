   ! --------------------------------------------------------------- !
   subroutine Base_Initialization(this, calculate_jacobian, provide_jacobian)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(inout) :: this
      logical, intent(in) :: calculate_jacobian, provide_jacobian
      ! ------------------------------------------------------------ !
      this%calculate_jacobian = calculate_jacobian
      this%provide_jacobian = provide_jacobian
      this%direct_variables = 0.0_dp
      this%direct_variables_mask = .False.
   end subroutine Base_Initialization
