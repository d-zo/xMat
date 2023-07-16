   ! --------------------------------------------------------------- !
   subroutine Base_Initialization(this, calculateJacobian, provideJacobian)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(inout) :: this
      logical, intent(in) :: calculateJacobian, provideJacobian
      ! ------------------------------------------------------------ !
      this%calculateJacobian = calculateJacobian
      this%provideJacobian = provideJacobian
      this%direct_variables = 0.0_dp
      this%direct_variables_mask = .False.
   end subroutine Base_Initialization
