   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Nonzero_Division, Norm, Trace, Squared, Deviatoric_Part, Dyadic_Product22, &
                                 Double_Contraction42, const_identity4d_sym, Pack_States, Unpack_States
      !
      class(Hypoplasticity_Wu92), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, N_mat, D_dir
      real(dp), dimension(__tensor__) :: L_mat
      real(dp) :: trT, cur_voidratio, dot_voidratio
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(__tensor__) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! change of void ratio `\dot{e} = \left(1+e_i\right)\tr{(\mathbf{D})}`

      trT = Trace(cur_stress)
      if (this%Is_Valid_Stress_State(stress=cur_stress) .and. (abs(trT) > setting_epsilon)) then
         L_mat = this%param_C1*trT*const_identity4d_sym &            ! `\mathcal{L} = C_1 \tr{(\mathbf{T})}\mathcal{I}^{id} + \frac{C_2}{\tr{(\mathbf{T})}}\mathbf{T}\otimes\mathbf{T}`
               + this%param_C2/trT*Dyadic_Product22(cur_stress, cur_stress)
         N_mat = (this%param_C3*Squared(cur_stress) &                ! `\mathbf{N} = \frac{C_3\mathbf{T}^2 + C_4\mathbf{\hat{T}}^2}{\tr{(\mathbf{T})}}`
               + this%param_C4*Squared(Deviatoric_Part(cur_stress)))/trT
         dot_stress = Double_Contraction42(L_mat, dot_strain) &      ! `\overset{\circ}{\mathbf{T}} = \mathcal{L}:\mathbf{D} + \mathbf{N}||\mathbf{D}||`
                    + N_mat*Norm(dot_strain)

         if ((this%calculateJacobian) .and. (.not. setting_numerical_jacobian)) then
            D_dir = Nonzero_Division(val=dot_strain, fac=Norm(dot_strain))
            dot_jac_stress = L_mat + Dyadic_Product22(N_mat, D_dir)
         end if
      else
         dot_stress = 0.0_dp
      end if

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
