   ! --------------------------------------------------------------- !
   pure function Get_Direct_Variables_Mask(this) result(direct_var_mask)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(in) :: this
      logical, dimension(setting_num_statevariables) :: direct_var_mask
      ! ------------------------------------------------------------ !
      direct_var_mask = this%direct_variables_mask
   end function Get_Direct_Variables_Mask
