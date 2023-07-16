   ! --------------------------------------------------------------- !
   subroutine Calculate_Results(this, exportstress, exportstates, exportjacobian, exportdt)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Is_Nan, Write_Warning, setting_epsilon, &
                                  setting_restrict_initial_substep, setting_initial_substep_scale
      use Debug, only: Formatval
      use Math_Operations, only: Nonzero_Division, Norm, Pack_States, Unpack_States
      !
      class(Xmat), intent(inout) :: this
      real(dp), dimension(__matrix__), intent(out) :: exportstress
      real(dp), dimension(this%nstates), intent(out) :: exportstates
      real(dp), dimension(__tensor__), intent(out) :: exportjacobian
      real(dp), intent(out) :: exportdt
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: cur_state, new_state
      real(dp), dimension(setting_num_statevariables) :: statevariables, new_statevariables
      real(dp), dimension(__tensor__) :: jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables
      integer, dimension(3) :: solver_statistics
      real(dp) :: norm_dot_strain
      logical :: normal_integration, successful_normal_integration

      jac_stress = 0.0_dp
      jac_statevariables = 0.0_dp
      solver_statistics = 0

      norm_dot_strain = Norm(this%assigned_dot_strain)
      normal_integration = ((norm_dot_strain > setting_epsilon) .and. (this%assigned_dt > setting_epsilon))
      successful_normal_integration = .False.

      if (setting_restrict_initial_substep) then
         !exportdt = min(Nonzero_Division(val=setting_initial_substep_scale, fac=norm_dot_strain), &
         !   this%assigned_dt)
         exportdt = this%assigned_dt / max(1, &
            idnint(norm_dot_strain*this%assigned_dt/setting_initial_substep_scale))
      else
         exportdt = this%assigned_dt
      end if

      ! --- Packing all cur_states in a long vector
      statevariables = this%statevariables
      cur_state = Pack_States(stress=this%saved_stress, jac_stress=jac_stress, statevariables=statevariables, &
         jac_statevariables=jac_statevariables)

      if (normal_integration) then
         ! --- Do the actual integration of the constitutive routine
         call this%xmat_solver%Integrate(integ_obj=this%xmat_constitutive_model, cur_state=cur_state, &
            cur_time=this%assigned_time, dt=this%assigned_dt, new_state=new_state, mod_dt=exportdt, &
            successful_integration=successful_normal_integration)
         solver_statistics = this%xmat_solver%Get_Statistics()
      end if

      if (successful_normal_integration) then
         ! --- Unpacking of all states out of a long vector
         call Unpack_States(input_states=new_state, stress=exportstress, jac_stress=exportjacobian, &
            statevariables=new_statevariables, jac_statevariables=jac_statevariables)
      else
         ! If either the norm of the dot strain or the timeincrement is (almost) zero or the integration fails,
         ! make one call to Get_Dot_State() of the constitutive model to get the jacobian
         new_state = this%xmat_constitutive_model%Get_Dot_State(ref_dt=this%assigned_dt, &
            cur_time=this%assigned_time, cur_state=cur_state)
         ! If the time increment is greater than zero, perform a single explicit Euler step
         if (this%assigned_dt > setting_epsilon) then
            new_state = cur_state + new_state*this%assigned_dt
         end if

         ! --- Unpacking of all states out of a long vector
         call Unpack_States(input_states=new_state, stress=exportstress, jac_stress=exportjacobian, &
            statevariables=new_statevariables, jac_statevariables=jac_statevariables)
         exportstress = this%saved_stress
         new_statevariables = this%statevariables

         ! No recommendation to change the time increment (in case integration failed and the value was changed)
         exportdt = this%assigned_dt
      end if

      exportstates = new_statevariables(1:this%nstates)
      ! NOTE: Uncomment the following to add the solver statistics to exportstates. Make sure that there are enough
      ! unused state variables because they are currently overwritten without any checks
      !exportstates(this%nstates-size(solver_statistics)+1:this%nstates) = real(solver_statistics, dp)

      exportjacobian = Nonzero_Division(val=exportjacobian, fac=this%assigned_dt)

      if (any(Is_Nan(exportstress))) then
         call Write_Warning('Calculate_Results: NaN entries found in stresses. Program will probably fail shortly' &
            // char(10) // Formatval('T', exportstress))
      end if

      if (any(Is_Nan(exportjacobian))) then
         call Write_Warning('Calculate_Results: NaN entries found in jacobian. ' &
            // 'Convergence with implicit algorithm unlikely' // char(10) // Formatval('jacobian', exportjacobian))
      end if

      if (any(Is_Nan(exportstates))) then
         call Write_Warning('Calculate_Results: NaN entries found in states. Program will probably fail shortly' &
            // char(10) // Formatval('states', exportstates))
      end if
   end subroutine Calculate_Results
