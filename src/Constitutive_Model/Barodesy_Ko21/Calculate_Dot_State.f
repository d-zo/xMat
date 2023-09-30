   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: const_identity2d, Nonzero_Division, Norm, Trace, Deviatoric_Part, &
                                 Matrix_Exponential, Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Barodesy_Ko21), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, stress_dir, R_mat, R_dir
      real(dp), dimension(__matrix__) :: D_dir, D_dirold, D_dirdev
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(__tensor__) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: cur_voidratio, dot_voidratio, p_mean, norm_D, norm_T, delta, e_c, lambda_D, &
                  dot_e_c1, dot_e_c2, eps1, eps2, h_fac, f_fac, g_fac

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      ! Currently e_c is saves in two fields due to different integration of both components
      e_c = statevariables(2) + statevariables(3)
      D_dirold = Vec9_To_Mat(vec9=statevariables(4:12))

      norm_D = Norm(dot_strain)
      D_dir = Nonzero_Division(val=dot_strain, fac=norm_D)
      delta = Trace(D_dir)

      norm_T = Norm(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)
      p_mean = -Trace(cur_stress)/3.0_dp

      eps1 = this%param_k1 - this%param_k2*log(p_mean)               ! (1) of Kolymbas (2021): `\epsilon_1(p) = k_1 - k_2 \log{(p)}`
      eps2 = eps1*(1.0_dp + this%param_kappa2*delta)                 ! (2) of Kolymbas (2021): `\epsilon_2(p) = \epsilon_1(p)\cdot \left(1 + \kappa_2 \delta\right)`

      D_dirdev = Deviatoric_Part(D_dir)
      R_mat = -Matrix_Exponential(mat=D_dirdev * this%derived_c1) &  ! (11) of Kolymbas (2021): `\mathbf{R}(\mathbf{D}) = -\exp{\left(c_1 \hat{\mathbf{D}}^0\right)} + c_2 \delta \mathbf{1}`
            + this%param_c2*delta*const_identity2d
      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))

      h_fac = this%param_c4 * norm_T**this%param_c5                  ! (12) of Kolymbas (2021): `h = c_4 ||\mathbf{T}||^{c_5}`
      f_fac = e_c + this%param_c3*delta                              ! `f = e_c + c_3 \delta`
      g_fac = -cur_voidratio + this%param_c3*delta                   ! `g = -e + c_3 \delta`

      dot_stress = h_fac * (f_fac*R_dir + g_fac*stress_dir) * norm_D ! (10) of Kolymbas (2021): `\overset{\circ}{\mathbf{T}} = h\cdot\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)

      ! The integration of dot_e_c relies on two different parts as hinted in Kolymbas (2021)
      dot_e_c1 = this%param_kappa1*(eps2 - e_c)                      ! First part of (3) of Kolymbas (2021): `\dot{e}_{c,1} = \kappa_1 \cdot \left(\epsilon_2(p) - e_c\right)\cdot \dot{\varepsilon}`
      !                                                              ! which will be integrated normally
      lambda_D = Norm(D_dir - D_dirold)
      dot_e_c2 = lambda_D * 0.5_dp * (this%param_e_c0 - e_c)         ! Second part of (3) of Kolymbas (2021): `\dot{e}_{c,2} =  \left(e_{c0} - e_c\right) \cdot \dot{\lambda_D} / 2`
      !                                                              ! which has to be added to the integrated `\dot{e}_{c,1}`

      this%direct_variables(3) = dot_e_c2
      this%direct_variables(4:12) = Mat_To_Vec9(mat=D_dir)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_statevariables(2) = dot_e_c1
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
