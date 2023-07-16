   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_min_youngs_modulus, Write_Error_And_Exit
      use Math_Operations, only: const_identity2d, const_identity4d_sym, const_identity4d_tr, &
                                 Nonzero_Division, Trace, Dimensionless, Deviatoric_Part, Norm, Inverse_Tensor, &
                                 Tensor_Partialtrace, Double_Contraction42, Double_Contraction22, Dyadic_Product22, &
                                 Double_Contraction44, Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Viscohypoplasticity_Ni03), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, T_dev, T_dless, T_dless_dev, B_mat, B_dir, D_vis, &
                                dpressure_preconsol_ds, N_mat, igran_strain, dot_igran_strain, N_hat
      real(dp), dimension(__tensor__) :: L_hat, L_mat, L_hat_inv, A_term, B_term, C_mat, C_inv, &
                                   K_firstinv, K_first, jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: trT, cur_voidratio, dot_voidratio, factor_F, factor_F2, state_OCR, &
                  cur_param_young, p_mean, q_mean, crit_stress_ratio, eta, eta2, &
                  dtth, pressure_equiv, pressure_preconsol, dpressure_preconsol_dp, dpressure_preconsol_dq
      logical :: successful_inversion

      ! ATTENTION: This routine is not fully tested and the current implementation will probably not work as intended.
      !            Do not use it unless you have checked it thouroughly and improved it first!

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

      if (.not. this%Valid_State(cur_stress=cur_stress, cur_voidratio=cur_voidratio, &
         pressure_equiv=pressure_equiv)) then

         if (cur_param_young < setting_min_youngs_modulus) then
            cur_param_young = setting_min_youngs_modulus
         end if

         call this%Elasticity(youngs_modulus=cur_param_young, nu=this%param_nu, dot_strain=dot_strain, &
            dot_stress=dot_stress, jacobian=dot_jac_stress)

         dot_voidratio = 0.0_dp                                      ! Keep current void ratio
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
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! change of void ratio `\dot{e} = \left(1+e_i\right)\tr{(\mathbf{D})}`

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      T_dev = Deviatoric_Part(cur_stress)                            ! Deviatoric part of the stress tensor `\mathbf{T}^\mathrm{dev}`
      T_dless = Dimensionless(cur_stress)                            ! Dimensionless stress tensor `\hat{\mathbf{T}}`
      T_dless_dev = Deviatoric_Part(T_dless)                         ! Deviatoric part of dimensionless stress tensor `\hat{\mathbf{T}}^\mathrm{dev}`

      ! Compute tensor `\mathcal{L}` and matrix `\mathbf{N}`
      factor_F = this%Get_Factor_F(stress_dless_dev=T_dless_dev)     ! (2.66) of Niemunis (2003): `F=\sqrt{\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
      factor_F2 = factor_F**2

      L_hat = (factor_F2 + this%derived_b2)*const_identity4d_sym &   ! (2.63) of Niemunis (2003) with increased shear stiffness from p.167:
            + this%derived_a2*Dyadic_Product22(T_dless, T_dless) &   ! `\hat{\mathcal{L}}_\mathrm{incr} = \left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}`
            - this%derived_b2/3.0_dp*const_identity4d_tr

      L_hat_inv = 1.0_dp/factor_F2*(const_identity4d_sym &           ! (4.140) of Niemunis (2003): `\hat{\mathcal{L}}^{-1} = \frac{1}{F^2}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right)`
                - Dyadic_Product22(T_dless, T_dless)/(factor_F2/this%derived_a2 &
                + Double_Contraction22(T_dless, T_dless)))
      N_hat = factor_F*this%derived_a*(T_dless + T_dless_dev)
      B_mat = Double_Contraction42(L_hat_inv, N_hat)                 ! (4.141) of Niemunis (2003): `\mathbf{B} = \hat{\mathcal{L}}^{-1}:\hat{\mathbf{N}} = \frac{a}{F}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right):\left(\hat{\mathbf{T}}+\hat{\mathbf{T}}^\mathrm{dev}\right)`

      ! Alternatively used
      ! f_1 = factor_F2 + this%derived_a2*Double_Contraction22(T_dless, T_dless)
      ! f_2 = -this%derived_a2*Double_Contraction22(T_dless, T_dless+T_dless_dev)
      ! B_mat = f_1*(T_dless + T_dless_dev) + f_2 * T_dless

      B_dir = Nonzero_Division(val=B_mat, fac=Norm(B_mat))

      p_mean = -trT/3.0_dp
      q_mean = sqrt(1.5_dp*Double_Contraction22(T_dev, T_dev))

      crit_stress_ratio = factor_F * this%derived_pq_incl            ! p.126 of Niemunis (2003): `M(\mathbf{T}) = F(\mathbf{T}) \frac{6\sin{(\varphi_c)}}{3-\sin{(\varphi_c)}}`
      eta = q_mean/(crit_stress_ratio*p_mean)                        ! `\bar{\eta} = \frac{q}{Mp}`
      eta2 = eta**2                                                  ! Within (4.80) of Niemunis (2003): `\bar{\eta}^2 =\left(\frac{q}{Mp}\right)^2`

      pressure_preconsol = this%Get_Pressure_Preconsol( &            ! `p_e^+`
                           overcrit=eta2, trT=trT)
      state_OCR = pressure_equiv / pressure_preconsol                ! (4.78) of Niemunis (2003): `OCR = \frac{p_e}{p_e^+}`
      L_mat = -this%derived_beta_b*trT*L_hat                         ! (4.64) and (4.66) of Niemunis (2003): `\mathcal{L} = f_b\hat{\mathcal{L}} = -\beta_b\tr{(\mathbf{T})}\hat{\mathcal{L}}`

      D_vis = -this%param_D_r*B_dir &                                ! (4.74) of Niemunis (2003): `\mathbf{D}^\mathrm{vis} = -D_r \vec{\mathbf{B}}\left(\frac{1}{OCR}\right)^{(1/I_v)}`
            * (1.0_dp/state_OCR)**(1.0_dp/this%param_I_v)

      if (eta > 1.0_dp) then                                         ! State above critical state surface `||\mathbf{B}|| = 1`
         pressure_preconsol = p_mean*(1.0_dp + eta2)/2.0_dp &        ! (4.118) of Niemunis (2003): `p_e^{+\mathrm{new}} = p\left(1+\bar{\eta}^2\right)\frac{1+\beta_B}{2}`
                            * (1.0_dp + this%param_beta_b)
         dpressure_preconsol_dp = (1.0_dp - eta2)/2.0_dp &           ! Modification of (4.119) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial p} = \frac{p_e^{+\mathrm{new}}}{p}`
                                * (1.0_dp + this%param_beta_b)
         dpressure_preconsol_dq = eta/crit_stress_ratio &            ! Modification of (4.120) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial q} = \frac{\bar{\eta}(\beta_B + 1)}{M}`
                                * (1.0_dp + this%param_beta_b)
      else                                                           ! State below critical state surface `||\mathbf{B}|| = 1`
         pressure_preconsol = p_mean*(1.0_dp - this%param_beta_b &   ! (4.117) of Niemunis (2003): `p_e^{+\mathrm{new}} = \frac{p}{\beta_B - 1}\left[\beta_B\sqrt{1+\bar{\eta}^2(\beta_B^2 - 1)} - 1\right]`
                            * sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp)))/(1.0_dp - this%param_beta_b)
         dpressure_preconsol_dp = pressure_preconsol/p_mean &        ! (4.119) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial p} = \frac{p_e^{+\mathrm{new}}}{p} - \frac{\beta_B \bar{\eta}^2(\beta_B + 1)}{\sqrt{1 + \bar{\eta}^2(\beta_B^2-1)}}`
                                - this%param_beta_b*eta2 * (this%param_beta_b + 1.0_dp) &
                                / sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp))
         dpressure_preconsol_dq = this%param_beta_b*eta &            ! (4.120) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial q} = \frac{\beta_B \bar{\eta}(\beta_B + 1)}{M\sqrt{1 + \bar{\eta}^2(\beta_B^2-1)}}`
                                / crit_stress_ratio*(this%param_beta_b + 1.0_dp) &
                                / sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp))
      end if

      ! Calculate `\frac{\partial p_e^{+\mathrm{new}}}{\partial \mathbf{T}}` as in (4.133) of Niemunis (2003)
      dpressure_preconsol_ds = -const_identity2d*dpressure_preconsol_dp/3.0_dp
      if (q_mean > 0.0001_dp) then
         dpressure_preconsol_ds = dpressure_preconsol_ds + T_dev*1.5_dp/q_mean*dpressure_preconsol_dq
      end if

      ! NOTE: Check dependence of dtth on ref_dt and compare it with other solutions. Although using dtth=0 would make
      !       the following procedure easier (L_mat and D_vis will stay), it might be critical to have a time dependence
      dtth = 0.5*ref_dt

      A_term = Dyadic_Product22(D_vis, dpressure_preconsol_ds)*dtth &! Compare to (4.132) of Niemunis (2003):
             * state_OCR/(this%param_I_v*pressure_equiv)             ! `2\mathcal{A} \approx \mathbf{D}^\mathrm{vis} \frac{OCR}{I_v \cdot p_e}\frac{\partial p_e^{+\mathrm{new}}}{\partial \mathbf{T}}`
      B_term = Dyadic_Product22(D_vis, const_identity2d)*dtth &      ! Compare to (4.134) of Niemunis (2003):
             * (1.0_dp + cur_voidratio) &                            ! `2\mathcal{B} = \frac{1+e^t}{I_v \lambda}\mathbf{D}_t^\mathrm{vis} \mathbf{1}`
             / (this%param_I_v*this%param_lambda)
      C_inv = const_identity4d_sym - B_term                          ! Inverse of (4.130) of Niemunis (2003): `\mathcal{C} = \left[\mathcal{M}^t-\mathcal{L}^t:\mathcal{B}\Delta t\right]^{-1}:\mathcal{L}^t`

      call Inverse_Tensor(tensor=C_inv, inv_tensor=C_mat, successful_inversion=successful_inversion)
      if (.not. successful_inversion) then
         call Write_Error_And_Exit('Viscohypoplasticity: Inversion of C_inv failed')
      end if

      K_firstinv = Double_Contraction44(L_mat, A_term) &             ! Inverse of first part of (4.129) of Niemunis (2003): `\mathcal{I} + \mathcal{L}^t:\mathcal{A}\Delta t`
                 + const_identity4d_sym

      call Inverse_Tensor(tensor=K_firstinv, inv_tensor=K_first, successful_inversion=successful_inversion)
      if (.not. successful_inversion) then
         call Write_Error_And_Exit('Viscohypoplasticity: Inversion of K_firstinv failed')
      end if

      L_mat = Double_Contraction44( &                                ! Similar to (4.129) of Niemunis (2003):
              Double_Contraction44(K_first, L_mat), C_inv)           ! `\mathcal{K} = \left[\mathcal{I} + \mathcal{L}^t:\mathcal{A}\Delta t\right]^{-1}:\left[\mathcal{M}^t - \mathcal{L}^t:\mathcal{B}\Delta t\right]`

      D_vis = Double_Contraction42(C_mat, D_vis)                     ! Apply transformation on `\mathbf{D}^\mathrm{vis}` as as below (4.131) of Niemunis (2003)

      N_mat = 0.0_dp
      call this%Intergranular_Strain(L_mat=L_mat, fdN_mat=N_mat, D_vis=D_vis, hypoplastic=.False., &
         igran_strain=igran_strain, R_max=this%param_R_max, m_T=this%param_m_T, m_R=this%param_m_R, &
         beta_R=this%param_beta_R, chi=this%param_chi, dt=ref_dt, dot_strain=dot_strain, &
         dot_stress=dot_stress, dot_igran_strain=dot_igran_strain, jacobian=dot_jac_stress)

      ! NOTE: Estimation of current stiffness for replacement model should be checked for correctness
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
