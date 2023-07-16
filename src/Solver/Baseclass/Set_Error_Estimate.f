   ! --------------------------------------------------------------- !
   subroutine Set_Error_Estimate(this, estimation)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      real(dp), dimension(setting_max_internal_states), intent(in) :: estimation
      ! ------------------------------------------------------------ !
      ! Mask all values (set to zero) which should not contribute to error estimation
      where (this%mask_ignored_indices)
         this%current_error = 0.0_dp
      elsewhere
         this%current_error = estimation
      end where
      this%current_error_assigned = .True.
   end subroutine Set_Error_Estimate
