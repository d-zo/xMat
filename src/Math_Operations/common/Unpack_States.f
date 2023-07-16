   ! --------------------------------------------------------------- !
   pure subroutine Unpack_States(input_states, stress, jac_stress, statevariables, jac_statevariables)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, setting_max_internal_states
      !
      real(dp), dimension(setting_max_internal_states), intent(in) :: input_states
      real(dp), dimension(__matrix__), intent(out) :: stress
      real(dp), dimension(__tensor__), intent(out) :: jac_stress
      real(dp), dimension(setting_num_statevariables), intent(out) :: statevariables
      real(dp), dimension(setting_num_statevariables, __matrix__), intent(out) :: jac_statevariables
      ! ------------------------------------------------------------ !
      stress = reshape(input_states(1:__nelmat__), [__matrix__])
      jac_stress = reshape(input_states(1+__nelmat__:(__nelmat__+1)*__nelmat__), [__tensor__])
      statevariables = input_states(1+(__nelmat__+1)*__nelmat__:(__nelmat__+1)*__nelmat__+setting_num_statevariables)
      jac_statevariables = reshape( &
         input_states(1+(__nelmat__+1)*__nelmat__+setting_num_statevariables:(__nelmat__+1)*__nelmat__+(__nelmat__+1)*setting_num_statevariables), &
         [setting_num_statevariables, __matrix__])
   end subroutine Unpack_States
