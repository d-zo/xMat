   ! --------------------------------------------------------------- !
   subroutine Integrate(this, integ_obj, cur_state, cur_time, dt, new_state, mod_dt, successful_integration)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Warning, &
                                  setting_stepsize_scaling_safety, setting_max_rel_error, &
                                  setting_max_integration_steps, setting_max_integration_refinements
      use Debug, only: Formatval
      !
      class(Solver), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states), intent(out) :: new_state
      real(dp), intent(inout) :: mod_dt
      logical, intent(out) :: successful_integration
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k1, inp_state, cum_error, state_scaled_err
      real(dp), dimension(setting_num_statevariables) :: tmp_direct_variables, direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask
      real(dp) :: steptime, stepsize, temp_error, temp_stepsize, next_stepsize
      integer :: idx_integloop, idx_convloop

      successful_integration = .False.
      this%last_integ_success = .False.
      this%last_num_steps_accepted = 0
      this%last_num_steps_rejected = 0
      cum_error = 0.0_dp
      inp_state = cur_state

      direct_variables = 0.0_dp

      idx_integloop = 0
      idx_convloop = 0

      steptime = 0.0_dp
      stepsize = mod_dt
      if (stepsize > dt) then
         ! Initial mod_dt should not be larger than overall dt
         stepsize = dt
      end if

      if ((this%stepsize_fixed) .and. (stepsize < dt/setting_max_integration_steps)) then
         call Write_Warning('Integrate: Fixed step integrator needs more steps (' // &
            Formatval('', int(dt/stepsize)) // ') than allowed ' &
            // 'in setting_max_integration_steps (' // Formatval('', setting_max_integration_steps) // ')')
         new_state = cur_state
         return
      end if

      integration_loop: &
      do
         idx_integloop = idx_integloop + 1

         if (this%next_dot_state_assigned) then
            dot_state_k1 = this%next_dot_state
            this%next_dot_state_assigned = .False.
            this%next_dot_state = 0.0_dp
         else
            dot_state_k1 = integ_obj%Get_Dot_State( &                ! `k_1 = f(t_n, y_n)`
               ref_dt=stepsize, cur_time=cur_time+steptime, cur_state=inp_state)
         end if
         call integ_obj%Get_Direct_Variables(direct_variables=tmp_direct_variables)

         ! Scaled error estimation from Press et al. (1997): Numerical Recipes in Fortran, p. 711/715
         state_scaled_err = abs(inp_state) + abs(stepsize*dot_state_k1) + setting_epsilon
         next_stepsize = stepsize

         idx_convloop = 0                                            ! Reset before entering convergence_loop
         convergence_loop: &
         do
            idx_convloop = idx_convloop + 1
            stepsize = next_stepsize

            new_state = this%Integration_Algorithm(integ_obj=integ_obj, cur_state=inp_state, &
               dot_state_k1=dot_state_k1, cur_time=cur_time + steptime, dt=stepsize)
            ! Check if stepsize is not fixed, error information is present and error is too big to reject step
            if (this%stepsize_fixed) then
               if (this%current_error_assigned) then
                  cum_error = cum_error + this%current_error
               end if

               this%last_num_steps_accepted = this%last_num_steps_accepted + 1
               exit convergence_loop
            end if

            if (.not. this%current_error_assigned) then
               this%last_num_steps_accepted = this%last_num_steps_accepted + 1
               exit convergence_loop
            end if

            temp_error = maxval(abs(this%current_error)/state_scaled_err)/setting_max_rel_error
            if (temp_error <= 1.0_dp) then
               ! Step is accepted and scaling is performed on next step size
               cum_error = cum_error + this%current_error
               this%last_num_steps_accepted = this%last_num_steps_accepted + 1

               next_stepsize = setting_stepsize_scaling_safety*stepsize*(temp_error**this%exp_stepgrow)
               next_stepsize = min(next_stepsize, this%max_stepgrow*stepsize)
               exit convergence_loop
            else
               ! Step is rejected. Step size for this step will be scaled
               this%last_num_steps_rejected = this%last_num_steps_rejected + 1

               temp_stepsize = setting_stepsize_scaling_safety*stepsize*(temp_error**this%exp_stepshrink)
               next_stepsize = max(temp_stepsize, this%min_stepshrink*stepsize)
            end if

            if (idx_convloop >= setting_max_integration_refinements) then
               call Write_Warning('Integrate: Too many iterations in convergence_loop')
               return
            end if

            if (next_stepsize <= setting_epsilon) then
               call Write_Warning('Integrate: Stepsize got reduced below minimum value')
               return
            end if
         end do convergence_loop

         direct_variables = direct_variables + tmp_direct_variables*stepsize/dt

         inp_state = new_state

         steptime = steptime + stepsize
         stepsize = next_stepsize
         mod_dt = stepsize

         if (steptime >= dt - setting_epsilon) then
            exit integration_loop
         end if

         if (steptime + stepsize > dt) then
            stepsize = dt - steptime
         end if

         if (idx_integloop >= setting_max_integration_steps) then
            call Write_Warning('Integrate: Too many iterations in integration_loop')
            return
         end if
      end do integration_loop

      if (this%current_error_assigned) then
         this%current_error = cum_error
      end if
      successful_integration = .True.
      this%last_integ_success = .True.

      ! Assign direct variables to their positions in new_state
      direct_variables_mask = integ_obj%Get_Direct_Variables_Mask()
      where (direct_variables_mask)
         new_state((__nelmat__+1)*__nelmat__+1:(__nelmat__+1)*__nelmat__+setting_num_statevariables) = direct_variables
      end where
   end subroutine Integrate
