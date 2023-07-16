   ! --------------------------------------------------------------- !
   subroutine Consistent_Jacobian(this, ref_dt, cur_time, cur_state, dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian_disturbance
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Perturbate, &
                                 Get_Matrix_From_Tensorpart, Get_Elements_From_Matrixlist, &
                                 Set_Matrix_In_Tensorpart, Set_Elements_In_Matrixlist
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states), intent(inout) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress, cur_stress_pert, dot_stress_pert
      real(dp), dimension(__matrix__) :: dot_strain_pert, tmp_jac_mat
      real(dp), dimension(setting_num_statevariables) :: cur_statevariables, dot_statevariables, tmp_jac_vec, &
                                                         cur_statevariables_pert, dot_statevariables_pert
      real(dp), dimension(__tensor__) :: cur_jacobian, dot_jacobian, tmp_jacobian
      real(dp), dimension(setting_num_statevariables, __matrix__) :: cur_sta_jacobian, dot_sta_jacobian, tmp_sta_jacobian
      real(dp), dimension(setting_max_internal_states) :: tmp_state, dot_state_pert
      integer :: idx_loop, idx, jdx

      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=cur_jacobian, &
         statevariables=cur_statevariables, jac_statevariables=cur_sta_jacobian)
      call Unpack_States(input_states=dot_state, stress=dot_stress, jac_stress=dot_jacobian, &
         statevariables=dot_statevariables, jac_statevariables=dot_sta_jacobian)
      ! The values of dot_jacobian and dot_sta_jacobian are assumed to be initialized as zero in the constitutive model(s)

      do idx_loop = 1, __nelmat__
         call Ref_Index(ref_idx=idx_loop, idx=idx, jdx=jdx)

         tmp_jac_mat = Get_Matrix_From_Tensorpart(tens=cur_jacobian, idx=idx, jdx=jdx)
         tmp_jac_vec = Get_Elements_From_Matrixlist(matlist=cur_sta_jacobian, nel=setting_num_statevariables, &
            idx=idx, jdx=jdx)

         cur_stress_pert = cur_stress + setting_numerical_jacobian_disturbance * tmp_jac_mat
         cur_statevariables_pert = cur_statevariables + setting_numerical_jacobian_disturbance * tmp_jac_vec
         tmp_state = Pack_States(stress=cur_stress_pert, jac_stress=cur_jacobian, &
            statevariables=cur_statevariables_pert, jac_statevariables=cur_sta_jacobian)

         dot_strain_pert =  this%dot_strain
         call Perturbate(mat=dot_strain_pert, idx=idx, jdx=jdx, theta=setting_numerical_jacobian_disturbance)

         dot_state_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=tmp_state, dot_strain=dot_strain_pert)

         call Unpack_States(input_states=dot_state_pert, stress=dot_stress_pert, jac_stress=tmp_jacobian, &
            statevariables=dot_statevariables_pert, jac_statevariables=tmp_sta_jacobian)

         tmp_jac_mat = (dot_stress_pert - dot_stress) / setting_numerical_jacobian_disturbance
         tmp_jac_vec = (dot_statevariables_pert - dot_statevariables) / setting_numerical_jacobian_disturbance

         call Set_Matrix_In_Tensorpart(tens=dot_jacobian, idx=idx, jdx=jdx, mat=tmp_jac_mat)
         call Set_Elements_In_Matrixlist(matlist=dot_sta_jacobian, nel=setting_num_statevariables, &
            idx=idx, jdx=jdx, elemlist=tmp_jac_vec)
      end do

      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jacobian, statevariables=dot_statevariables, &
         jac_statevariables=dot_sta_jacobian)
   end subroutine Consistent_Jacobian
