   ! --------------------------------------------------------------- !
   pure function Pack_States(stress, jac_stress, statevariables, jac_statevariables) result(output_states)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, setting_max_internal_states
      !
      real(dp), dimension(__matrix__), intent(in) :: stress
      real(dp), dimension(__tensor__), intent(in) :: jac_stress
      real(dp), dimension(setting_num_statevariables), intent(in) :: statevariables
      real(dp), dimension(setting_num_statevariables, __matrix__), intent(in) :: jac_statevariables
      real(dp), dimension(setting_max_internal_states) :: output_states
      ! ------------------------------------------------------------ !
      output_states(1:__nelmat__) = reshape(stress, [__nelmat__])
      output_states(1+__nelmat__:(__nelmat__+1)*__nelmat__) = reshape(jac_stress, [__nelmat__*__nelmat__])
      output_states(1+(__nelmat__+1)*__nelmat__:(__nelmat__+1)*__nelmat__+setting_num_statevariables) = statevariables
      output_states(1+(__nelmat__+1)*__nelmat__+setting_num_statevariables:(__nelmat__+1)*__nelmat__+(__nelmat__+1)*setting_num_statevariables) = &
         reshape(jac_statevariables, [setting_num_statevariables*__nelmat__])
   end function Pack_States
