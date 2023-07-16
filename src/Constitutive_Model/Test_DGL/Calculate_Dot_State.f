   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Set_Element_In_Tensor, Double_Contraction42
      !
      class(Test_DGL), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(__matrix__), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(__matrix__) :: cur_stress, dot_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(__tensor__) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, __matrix__) :: jac_statevariables, dot_jac_statevariables
      !
      real(dp) :: assign_val
      integer :: ref_idx, ref_jdx, idx, jdx, kdx, ldx
      real(dp) :: mod_idx, mod_jdx
      logical :: symmetric

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)

      ! --- Calculate dot_stress (response of test DGL) based on the selected test case
      if (this%param_select == 1) then
         dot_stress = 10.0_dp * cur_stress
         ! Solution for this DGL is:   stress_0 * exp(10*t)
      else if (this%param_select == 2) then
         dot_stress = 0.8_dp*(exp(2.5_dp)-exp(0.5_dp))*exp(-0.25_dp*cur_stress)
         ! Solution for this DGL is:   4 * log(0.2*t*(exp(2.5)-exp(0.5)) + exp(0.5))
      else if (this%param_select == 3) then
         dot_stress = log(9.0_dp)/5.0_dp*(cur_stress - 1.0_dp)
         ! Solution for this DGL is:   9**(0.2*t) + 1
      else
         symmetric = .True.
         do ref_idx = 1, __nelmat__
            call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
            do ref_jdx = 1, __nelmat__
               call Ref_Index(ref_idx=ref_jdx, idx=kdx, jdx=ldx)
               assign_val = 10.0_dp*(ref_idx-1)+ref_jdx
               if (__nelmat__ == 6) then
                  mod_idx = real(ref_idx, dp)
                  if (ref_idx > 3) then
                     mod_idx = 2 + (mod_idx - 3)*2
                  end if
                  mod_jdx = real(ref_jdx, dp)
                  if (ref_jdx > 3) then
                     mod_jdx = 2 + (mod_jdx - 3)*2
                  end if
                  assign_val = 10.0_dp*(mod_idx-1)+mod_jdx
               else if (symmetric) then
                  mod_idx = real(ref_idx, dp)
                  if (ref_idx > 3) then
                     mod_idx = 2 + floor((mod_idx - 2)/2)*2
                  end if
                  mod_jdx = real(ref_jdx, dp)
                  if (ref_jdx > 3) then
                     mod_jdx = 2 + floor((mod_jdx - 2)/2)*2
                  end if
                  assign_val = 10.0_dp*(mod_idx-1)+mod_jdx
               end if
               call Set_Element_In_Tensor(tens=dot_jac_stress, idx=idx, jdx=jdx, kdx=kdx, ldx=ldx, val=assign_val)
            end do
            dot_stress = Double_Contraction42(dot_jac_stress, dot_strain)
         end do
      end if

      ! --- Packing all dot_states in a long vector
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
