   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_min_youngs_modulus, setting_default_nu
      use Math_Operations, only: Nonzero_Division, Norm, Trace, Matrix_Exponential, Pack_States, Unpack_States
      !
      class(Barodesy_Ko15), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, R_mat, R_dir, D_dir, stress_dir
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(__tensor__) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: cur_voidratio, dot_voidratio, norm_D, norm_T, dilatancy, B_fak, e_c, &
                  h_fac, f_fac, g_fac, cur_param_young

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      cur_param_young = statevariables(2)

      ! NOTE: Valid_State has to return .False. if cur_stress is almost zero
      if (.not. this%Valid_State(cur_stress=cur_stress)) then
         if (cur_param_young < setting_min_youngs_modulus) then
            cur_param_young = setting_min_youngs_modulus
         end if

         call this%Elasticity(youngs_modulus=cur_param_young, nu=setting_default_nu, dot_strain=dot_strain, &
            dot_stress=dot_stress, jacobian=dot_jac_stress)

         this%direct_variables(2) = cur_param_young

         ! --- Packing all dot_states in a long vector
         dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
            jac_statevariables=dot_jac_statevariables)
         return
      end if

      norm_D = Norm(dot_strain)
      D_dir = Nonzero_Division(val=dot_strain, fac=norm_D)
      dilatancy = Trace(D_dir)

      norm_T = Norm(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)

      R_mat = -Matrix_Exponential(mat=this%derived_c1*D_dir &        ! (9) of Kolymbas (2015): `\mathbf{R} = -e^{c_1\mathbf{D}^0 e^{(c_2 \delta)}}`
            * exp(this%param_c2*dilatancy))

      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))           ! `\mathbf{R}^0 = \frac{\mathbf{R}}{||\mathbf{R}||}`

      B_fak = (this%param_e_c0 - this%param_e_min) &                 ! (27) of Kolymbas (2015): `B = \frac{e_{c0} - e_\text{min}}{e_{c0} + 1} \left(\frac{c_4 + c_5\||\mathbf{T}||}{c_4}\right)^{-(1 + e_\text{min})/c_5}`
            / (this%param_e_c0 + 1.0_dp) &
            * ((this%param_c4 + this%param_c5*norm_T)/this%param_c4)**(-(1.0_dp + this%param_e_min)/this%param_c5)
      e_c = (this%param_e_min + B_fak)/(1.0_dp - B_fak)              ! (26) of Kolymbas (2015): `e_c = \frac{e_\text{min}+B}{1-B}`

      f_fac = dilatancy + this%param_c3*e_c                          ! (29) of Kolymbas (2015): `f = \delta + c_3 e_c`
      g_fac = -this%param_c3*cur_voidratio                           ! (30) of Kolymbas (2015): `g = -c_3 e`
      h_fac = -(this%param_c4 + this%param_c5*norm_T) &              ! (23) of Kolymbas (2015): `h = \frac{c_4 + c_5||\mathbf{T}||}{e - e_\text{min}}`
            / (cur_voidratio - this%param_e_min)

      dot_stress = h_fac*(f_fac*R_dir + g_fac*stress_dir)*norm_D     ! (12) of Kolymbas (2015): `\overset{\circ}{\mathbf{T}} = h\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
