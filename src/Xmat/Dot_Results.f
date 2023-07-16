   ! --------------------------------------------------------------- !
   subroutine Dot_Results(this, dotstress, dotstates, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Is_Nan, Write_Warning
      use Debug, only: Formatval
      use Math_Operations, only: Pack_States, Unpack_States
      !
      class(Xmat), intent(inout) :: this
      real(dp), dimension(__matrix__), intent(out) :: dotstress
      real(dp), dimension(this%nstates), intent(out) :: dotstates
      real(dp), dimension(__tensor__), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: cur_state, dot_state
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables, direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask
      real(dp), dimension(__tensor__) :: jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables

      jac_stress = 0.0_dp
      jac_statevariables = 0.0_dp

      ! --- Packing all cur_states in a long vector
      statevariables = this%statevariables
      cur_state = Pack_States(stress=this%saved_stress, jac_stress=jac_stress, statevariables=statevariables, &
         jac_statevariables=jac_statevariables)

      ! --- Calculate the response of the constitutive model
      dot_state = this%xmat_constitutive_model%Get_Dot_State(ref_dt=this%assigned_dt, cur_time=this%assigned_time, &
         cur_state=cur_state)
      call this%xmat_constitutive_model%Get_Direct_Variables(direct_variables=direct_variables)

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=dot_state, stress=dotstress, jac_stress=jacobian, &
         statevariables=dot_statevariables, jac_statevariables=jac_statevariables)
      direct_variables_mask = this%xmat_constitutive_model%Get_Direct_Variables_Mask()
      where (direct_variables_mask)
         dot_statevariables =  direct_variables
      end where
      dotstates = dot_statevariables(1:this%nstates)

      if (any(Is_Nan(dotstress))) then
         call Write_Warning('Dot_Results: NaN entries found in stresses.' // char(10) // Formatval('T', dotstress))
      end if

      if (any(Is_Nan(jacobian))) then
         call Write_Warning('Dot_Results: NaN entries found in jacobian. ' &
            // 'Convergence with implicit algorithm unlikely' // char(10) // Formatval('jacobian', jacobian))
      end if

      if (any(Is_Nan(dotstates))) then
         call Write_Warning('Dot_Results: NaN entries found in states.' // char(10) // Formatval('states', dotstates))
      end if
   end subroutine Dot_Results
