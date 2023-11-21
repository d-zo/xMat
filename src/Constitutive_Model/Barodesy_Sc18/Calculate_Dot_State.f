   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_min_youngs_modulus, setting_default_nu
      use Math_Operations, only: const_root3, Nonzero_Division, Norm, Trace, Matrix_Exponential, &
                                 Pack_States, Unpack_States
      !
      class(Barodesy_Sc18), intent(inout) :: this
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
      real(dp) :: cur_voidratio, dot_voidratio, alpha, norm_D, norm_T, dilatancy, trT, e_c, &
                  h_fac, b_interp, f_fac, g_fac, cur_param_young

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
      dilatancy = Trace(D_dir)                                       ! (1.7) of Schranz (2018): `\delta = \frac{\tr{(\mathbf{D})}}{||\mathbf{D}||}`

      norm_T = Norm(cur_stress)
      trT = Trace(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)

      alpha = -30.0_dp + this%derived_c3 &                           ! `\alpha = -30 + c_3\frac{|\delta - \sqrt{2}|^{c_2}}{\left(1+|\delta-\sqrt{2}|\right)^{c_1}}`
            * abs(dilatancy - sqrt(2.0_dp))**this%derived_c2/(1.0_dp + abs(dilatancy - sqrt(2.0_dp)))**this%derived_c1
      R_mat = -Matrix_Exponential(mat=D_dir*alpha)                   ! `\mathbf{R} = -e^{\alpha\mathbf{D}^0}`

      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))           ! `\mathbf{R}^0 = \frac{\mathbf{R}}{||\mathbf{R}||}`

      e_c = (1.0_dp + this%param_e_c0) &                             ! Modifed version of Appendix C of Schranz (2018):
          * exp(-(-trT/3.0_dp)**(1.0_dp - this%param_xi) &           ! `e_c = \left(1+e_{c0}\right)e^{\left(-\frac{p^{1-\xi}}{K_r\left(1-\xi\right)}\right)} - 1`
          / (this%param_K_r*(1.0_dp - this%param_xi))) - 1.0_dp

      h_fac = this%derived_c4*norm_T**this%param_xi                  ! (4.32) of Schranz (2018): `h = c_4 ||\mathbf{T}||^\xi`

      b_interp = (this%derived_b_ext - this%derived_b_comp) &        ! (4.57) of Schranz (2018): `b = \left(b_\mathrm{ext}-b_\mathrm{comp}\right)\left(\frac{\delta+\sqrt{3}}{2\sqrt{3}}\right)^{c_7} + b_\mathrm{comp}`
               * ((dilatancy + const_root3)/(2.0_dp*const_root3))**this%param_c7 + this%derived_b_comp

      f_fac = this%derived_c8*b_interp*dilatancy + this%param_c6     ! (4.59) of Schranz (2018): `f = c_8 b \delta + c_6`
      g_fac = (1.0_dp - this%derived_c8)*b_interp*dilatancy &        ! (4.59) of Schranz (2018): `g = (1-c_8) b \delta + c_5\left[\left(\frac{1+e}{1+e_c}\right)^\zeta - 1\right] -c_6`
            + this%param_c5*(((1.0_dp + cur_voidratio)/(1.0_dp + e_c))**this%derived_zeta - 1.0_dp) - this%param_c6

      dot_stress = h_fac*(f_fac*R_dir + g_fac*stress_dir)*norm_D     ! (3.14) of Schranz (2018): `\overset{\circ}{\mathbf{T}} = h\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! As (4.33) of Schranz (2018)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
