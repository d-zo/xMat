   ! --------------------------------------------------------------- !
   subroutine Set_State_Mask(this, statemask)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      logical, dimension(setting_num_statevariables), intent(in) :: statemask
      ! ------------------------------------------------------------ !
      this%mask_ignored_indices((__nelmat__+1)*__nelmat__+1:(__nelmat__+1)*__nelmat__+setting_num_statevariables) = statemask
   end subroutine Set_State_Mask
