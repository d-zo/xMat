   ! --------------------------------------------------------------- !
   subroutine Initialize(this, name, exp_stepgrow, max_stepgrow, exp_stepshrink, min_stepshrink, stepsize_fixed)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), optional :: exp_stepgrow
      real(dp), intent(in), optional :: max_stepgrow
      real(dp), intent(in), optional :: exp_stepshrink
      real(dp), intent(in), optional :: min_stepshrink
      logical, intent(in), optional :: stepsize_fixed
      ! ------------------------------------------------------------ !
      this%name = name

      this%current_error_assigned = .False.
      this%current_error = 0.0_dp

      this%next_dot_state_assigned = .False.
      this%next_dot_state = 0.0_dp

      this%mask_ignored_indices = .False.

      this%last_integ_success = .False.
      this%last_num_steps_accepted = 0
      this%last_num_steps_rejected = 0

      if (present(exp_stepgrow)) then
         this%exp_stepgrow = exp_stepgrow
      else
         this%exp_stepgrow = -0.2_dp
      end if

      if (present(exp_stepshrink)) then
         this%exp_stepshrink = exp_stepshrink
      else
         this%exp_stepshrink = -0.25_dp
      end if

      if (present(max_stepgrow)) then
         this%max_stepgrow = max_stepgrow
      else
         this%max_stepgrow = 5.0_dp
      end if

      if (present(min_stepshrink)) then
         this%min_stepshrink = min_stepshrink
      else
         this%min_stepshrink = 0.1_dp
      end if

      if (present(stepsize_fixed)) then
         this%stepsize_fixed = stepsize_fixed
      else
         this%stepsize_fixed = .False.
      end if
   end subroutine Initialize
