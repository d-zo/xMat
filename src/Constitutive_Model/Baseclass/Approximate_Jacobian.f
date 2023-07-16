   ! --------------------------------------------------------------- !
   subroutine Approximate_Jacobian(this, ref_dt, cur_time, cur_state, dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian_disturbance
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Perturbate, Set_Matrix_In_Tensorpart
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states), intent(inout) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: tmp_stress, dot_stress, dot_stress1_pert, dot_stress2_pert
      real(dp), dimension(__matrix__) :: dot_strain_pert, tmp_jac_mat
      real(dp), dimension(setting_num_statevariables) :: tmp_statevariables, dot_statevariables
      real(dp), dimension(__tensor__) :: tmp_jacobian, dot_jacobian
      real(dp), dimension(setting_num_statevariables, __matrix__) :: tmp_sta_jacobian, dot_sta_jacobian
      real(dp), dimension(setting_max_internal_states) :: dot_state1_pert, dot_state2_pert
      integer :: idx_loop, idx, jdx

      call Unpack_States(input_states=dot_state, stress=dot_stress, jac_stress=dot_jacobian, &
         statevariables=dot_statevariables, jac_statevariables=dot_sta_jacobian)
      ! The values of dot_jacobian and dot_sta_jacobian are assumed to be initialized as zero in Calculate_Results()

      do idx_loop = 1, __nelmat__
         call Ref_Index(ref_idx=idx_loop, idx=idx, jdx=jdx)

         dot_strain_pert = 0.0_dp
         call Perturbate(mat=dot_strain_pert, idx=idx, jdx=jdx, theta=-setting_numerical_jacobian_disturbance)

         dot_state1_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=cur_state, dot_strain=dot_strain_pert)
         call Unpack_States(input_states=dot_state1_pert, stress=dot_stress1_pert, jac_stress=tmp_jacobian, &
            statevariables=tmp_statevariables, jac_statevariables=tmp_sta_jacobian)

         dot_state2_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=cur_state, dot_strain=-dot_strain_pert)
         call Unpack_States(input_states=dot_state2_pert, stress=dot_stress2_pert, jac_stress=tmp_jacobian, &
            statevariables=tmp_statevariables, jac_statevariables=tmp_sta_jacobian)

         tmp_jac_mat= (dot_stress2_pert - dot_stress1_pert) / (2.0_dp*setting_numerical_jacobian_disturbance)
         call Set_Matrix_In_Tensorpart(tens=dot_jacobian, idx=idx, jdx=jdx, mat=tmp_jac_mat)
      end do

      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jacobian, statevariables=dot_statevariables, &
         jac_statevariables=dot_sta_jacobian)
   end subroutine Approximate_Jacobian
