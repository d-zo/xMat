   ! --------------------------------------------------------------- !
   subroutine Import_State(this, nstates, statevariables, stress, dot_strain, totaltime, timeincrement)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, Write_Error_And_Exit
      use Debug, only: Formatval
      !
      class(Xmat), intent(inout) :: this
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: statevariables
      real(dp), dimension(__matrix__), intent(in) :: stress, dot_strain
      real(dp), intent(in) :: totaltime, timeincrement
      ! ------------------------------------------------------------ !
      if (nstates > setting_num_statevariables) then
         call Write_Error_And_Exit('Import_State: ' // Formatval('num_states', nstates) // ' requested but only ' // &
            Formatval('setting_num_statevariables', setting_num_statevariables) // &
            ' provided - increase setting_num_statevariables?')
      end if
      this%statevariables = 0.0_dp
      this%statevariables(1:nstates) = statevariables
      this%nstates = nstates
      this%assigned_dot_strain = dot_strain
      this%assigned_time = totaltime
      this%assigned_dt = timeincrement

      call this%xmat_constitutive_model%Set_Values(overall_dt=this%assigned_dt, &
         dot_strain=this%assigned_dot_strain)

      this%saved_stress = stress
   end subroutine Import_State
