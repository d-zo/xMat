   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, setting_hypo_increased_stiffness, setting_hypo_consistent_f_d, &
                                  Write_Error_And_Exit, setting_min_youngs_modulus
      use Math_Operations, only: const_root3, const_identity4d_sym, const_identity4d_tr, &
                                 Nonzero_Division, Norm, Trace, Dimensionless, Deviatoric_Part, Inverse_Tensor, &
                                 Tensor_Partialtrace, Dyadic_Product22, Double_Contraction22, Double_Contraction42, &
                                 Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Hypoplasticity_VW96), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, T_dless, T_dless_dev
      real(dp) :: cur_voidratio, dot_voidratio, trT, e_d, e_c, e_i, norm_D, cur_param_young
      real(dp) :: m_dir_voidratio, m_dir_stress, factor_F, factor_F2, f_LN, f_e, f_b, f_d_bar, f_d_sign, f_d, yout
      real(dp), dimension(__matrix__) :: fdN_mat, igran_strain, dot_igran_strain, D_vis
      real(dp), dimension(__tensor__) :: L_mat, L_inv, jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables
      logical :: successful_inversion

      dot_igran_strain = 0.0_dp

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      igran_strain = Vec9_To_Mat(vec9=statevariables(2:10))
      cur_param_young = statevariables(11)

      ! NOTE: Valid_State has to return .False. if either cur_voidratio or cur_stress are zero
      !       because they are used in the denominator of some hypoplastic expressions below.
      !       Also dot_voidratio will be initialised in Valid_State. It will be zero if the
      !       current void_ratio is already in valid bounds or non-zero otherwise.
      if (.not. this%Valid_State(cur_stress=cur_stress, cur_voidratio=cur_voidratio, &
         e_d=e_d, e_c=e_c, e_i=e_i, dt=ref_dt, dot_voidratio=dot_voidratio)) then

         if (cur_param_young < setting_min_youngs_modulus) then
            cur_param_young = setting_min_youngs_modulus
         end if

         call this%Elasticity(youngs_modulus=cur_param_young, nu=this%param_nu, dot_strain=dot_strain, &
            dot_stress=dot_stress, jacobian=dot_jac_stress)

         dot_igran_strain = Nonzero_Division(val=-igran_strain, &    ! Reset intergranular strain for replacement model
            fac=ref_dt)

         this%direct_variables(11) = cur_param_young

         ! --- Packing all dot_states in a long vector
         dot_statevariables(1) = dot_voidratio
         dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
         dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
            jac_statevariables=dot_jac_statevariables)
         return
      end if

      dot_voidratio = dot_voidratio + &                              ! change of void ratio `\dot{e} = \left(1+e\right)\tr{(\mathbf{D})}`
                      (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! additionally take into account if `e` was not in valid bounds and
                                                                     ! add the necessary change `\dot{e}_{corr}`

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      T_dless = Dimensionless(cur_stress)                            ! Dimensionless stress tensor `\hat{\mathbf{T}}`
      T_dless_dev = Deviatoric_Part(T_dless)                         ! Deviatoric part of dimensionless stress tensor `\hat{\mathbf{T}}^\mathrm{dev}`

      norm_D = Norm(dot_strain)                                      ! Norm of strain increments `||\mathbf{D}||`

      ! Compute scalar factors `f_e`, `f_b` (both barotropy) and `f_d` (pycnotropy)
      f_e = (e_c/cur_voidratio)**this%param_beta_H                   ! (2.70) of Niemunis (2003): `f_e(\tr{(\mathbf{T})},e) = \left(\frac{e_c}{e}\right)^\beta_H`
      f_b = (this%param_e_i0/this%param_e_c0)**this%param_beta_H &   ! (2.72) of Niemunis (2003): `f_b(\tr{(\mathbf{T})}) = \left(\frac{e_{i0}}{e_{c0}}\right)^\beta_H \frac{h_s}{n_H} \frac{1+e_i}{e_i} \left(\frac{-\tr{(\mathbf{T})}}{h_s}\right)^{1-n_H} \left[f_{b,\mathrm{lastterm}}\right]^{-1}`
          * this%param_h_s/this%param_n_H * (1.0_dp + e_i)/e_i &     ! Therefore `f_b(\tr{(\mathbf{T})}) = \left(\frac{e_{i0}}{e_{c0}}\right)^\beta_H \frac{h_s}{n_H} \frac{1+e_i}{e_i} \left(\frac{-\tr{(\mathbf{T})}}{h_s}\right)^{1-n_H} \left[3+a^2-a\sqrt{3}\left(\frac{e_{i0}-e_{d0}}{e_{c0}-e_{d0}}\right)^{\alpha_H}\right]^{-1}`
          * (-trT/this%param_h_s)**(1.0_dp - this%param_n_H) / this%derived_f_b_lastterm

      ! Either use standard definition of `f_d` as in (2.71) or enforce a consistent lower bound as in (4.223) of Niemunis (2003)
      consistent_f_d : &
      if ((setting_hypo_consistent_f_d) &                            ! The changes of the lower bound would only be applied for `e < e_c`
         .and. (cur_voidratio < e_c)) then

         m_dir_voidratio = 1.0_dp                                    ! Direction of void ratio `M_e^{(d)}` and
         m_dir_stress = -e_d/this%param_h_s * this%param_n_H &       ! direction of stresses `M_T^{(d)}` as in (4.211) of Niemunis (2003)
                      * (-trT/this%param_h_s)**(this%param_n_H-1.0_dp)
         f_d_bar = -(m_dir_voidratio*const_root3 &                   ! (4.221) of Niemunis (2003): `\bar{f}_d = -\frac{M_e^{(d)}\sqrt{3}(1+e)+M_T^{(d)}f_b f_e \frac{3}{\sqrt{3}}(3+a^2)}{M_T^{(d)}f_b f_e 3a}`
                 * (1.0_dp + cur_voidratio) + m_dir_stress*f_b*f_e*const_root3 &
                 * (3.0_dp + this%derived_a2)) / (m_dir_stress*f_b*f_e*3.0_dp*this%derived_a)
         f_d_sign = sign(1.0_dp, cur_voidratio - e_d) &              ! Factor from (4.222) of Niemunis (2003)
                  * (abs(cur_voidratio - e_d) &                      ! `f_{d,\mathrm{sign}} = \left(\frac{e-e_d}{e_c-e_d}\right)^{\alpha_H}` if `e > e_d` and `-\left(\frac{|e-e_d|}{e_c-e_d}\right)^{\alpha_H}` otherwise
                  / (e_c - e_d))**this%param_alpha_H
         f_d = f_d_sign + (1.0_dp - f_d_sign)**5 * f_d_bar           ! (4.223) of Niemunis (2003): `f_d = f_{d,\mathrm{sign}} + \left[1 - f_{d,\mathrm{sign}}\right]^z\bar{f}_d` with `z=5`
      else
         f_d = ((cur_voidratio - e_d) / (e_c - e_d)) &               ! (2.71) of Niemunis (2003): `f_d(\tr{(\mathbf{T})},e) = \left(\frac{e-e_d}{e_c-e_d}\right)^{\alpha_H} = r_e^{\alpha_H}`
             **this%param_alpha_H
      end if consistent_f_d

      ! Compute tensor `\mathcal{L}` and matrix `\mathbf{N}`
      factor_F = this%Get_Factor_F(stress_dless_dev=T_dless_dev)     ! (2.66) of Niemunis (2003): `F=\sqrt{\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
      factor_F2 = factor_F**2
      f_LN = f_b*f_e/Double_Contraction22(T_dless, T_dless)          ! `f_\mathrm{LN} = \frac{f_b f_e}{\hat{\mathbf{T}}:\hat{\mathbf{T}}}`

      L_mat = f_LN*(factor_F2*const_identity4d_sym &                 ! (2.63) of Niemunis (2003): `\mathcal{L} = f_\mathrm{LN}a^2 \cdot \left(\left(\frac{F}{a}\right)^2 \mathcal{I} + \hat{\mathbf{T}}\otimes\hat{\mathbf{T}}\right)`
            + this%derived_a2*Dyadic_Product22(T_dless, T_dless))
      fdN_mat = f_d*f_LN*factor_F*this%derived_a &                   ! Adjusted (2.64) of Niemunis (2003): `f_d\mathbf{N} = f_d f_\mathrm{LN}a^2 \cdot \left(\frac{F}{a}\right)\left(\hat{\mathbf{T}} + \hat{\mathbf{T}}^\mathrm{dev}\right)`
              * (T_dless + T_dless_dev)

      if (setting_hypo_increased_stiffness) then                     ! Increase shear stiffness
         L_mat = L_mat + f_LN*this%derived_b2*const_identity4d_sym & ! Modification of (2.63) of Niemunis (2003):
               - this%derived_b2/3.0_dp*const_identity4d_tr          ! `\mathcal{L}_\mathrm{incr} = f_\mathrm{LN} \cdot \left(\left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}\right)`
         L_inv = (const_identity4d_sym - Dyadic_Product22(T_dless, & ! Adaption of (3.33) of Niemunis (2003):
                 T_dless)/(factor_F2/this%derived_a2 &               ! `\mathcal{L}^{-1} = \frac{1}{f_\mathrm{LN}F^2}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right)`
               + Double_Contraction22(T_dless, T_dless)))/(f_LN*factor_F2)
         fdN_mat = Double_Contraction42(L_mat, &                     ! Forge `f_d\mathbf{N}_\mathrm{incr} = f_d\mathcal{L}_\mathrm{incr}:\left(\mathcal{L}^{-1}:\mathbf{N}\right)`
                   Double_Contraction42(L_inv, fdN_mat))             ! to preserve flow direction
      end if

      if (this%param_m_R <= 2.0_dp) then
         ! ------ Hypoplasticity without intergranular strain ------ !
         dot_stress = Double_Contraction42(L_mat, dot_strain) &      ! (2.61) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{L}:\mathbf{D} + f_d \mathbf{N}||\mathbf{D}||`
                    + fdN_mat*norm_D
         if ((this%calculateJacobian) .and. (.not. setting_numerical_jacobian)) then
            if ((setting_hypo_increased_stiffness) .and. (cur_time > setting_epsilon) &
               .and. (norm_D > setting_epsilon)) then

               call Inverse_Tensor(tensor=L_inv, inv_tensor=L_mat, & ! Inverse of `\mathcal{L}_\mathrm{incr} = f_\mathrm{LN} \cdot \left(\left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}\right)`
                  successful_inversion=successful_inversion)
               if (.not. successful_inversion) then
                  call Write_Error_And_Exit('Hypoplasticity: Inversion of L_incr failed')
               end if

               yout = min(0.95_dp, 0.95_dp/(Norm(Double_Contraction42(L_inv, fdN_mat)) + setting_epsilon))
               dot_jac_stress = L_mat + yout*Dyadic_Product22(fdN_mat, dot_strain/norm_D)
            else
               dot_jac_stress = L_mat
            end if
         end if
      else
         ! -------- Hypoplasticity with intergranular strain ------- !
         D_vis = 0.0_dp
         call this%Intergranular_Strain(L_mat=L_mat, fdN_mat=fdN_mat, D_vis=D_vis, hypoplastic=.True., &
            igran_strain=igran_strain, R_max=this%param_R_max, m_T=this%param_m_T, m_R=this%param_m_R, &
            beta_R=this%param_beta_R, chi=this%param_chi, dt=ref_dt, dot_strain=dot_strain, &
            dot_stress=dot_stress, dot_igran_strain=dot_igran_strain, jacobian=dot_jac_stress)
      end if

      ! Estimation of current stiffness for replacement model
      if ((this%calculateJacobian) .and. (.not. setting_numerical_jacobian)) then
         cur_param_young = Tensor_Partialtrace(dot_jac_stress)/3.0_dp
      else
         cur_param_young = Tensor_Partialtrace(L_mat)/3.0_dp
      end if

      this%direct_variables(11) = cur_param_young

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
