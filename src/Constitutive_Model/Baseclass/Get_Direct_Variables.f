   ! --------------------------------------------------------------- !
   subroutine Get_Direct_Variables(this, direct_variables)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(in) :: this
      real(dp), dimension(setting_num_statevariables), intent(out) :: direct_variables
      ! ------------------------------------------------------------ !
      direct_variables = this%direct_variables
   end subroutine Get_Direct_Variables
