   ! --------------------------------------------------------------- !
   function Get_Dot_State(this, ref_dt, cur_time, cur_state) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian, setting_epsilon
      use Math_Operations, only: Norm
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      dot_state = this%Calculate_Dot_State(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_strain=this%dot_strain)

      if ((this%calculate_jacobian) .and. ((setting_numerical_jacobian) .or. (.not. this%provide_jacobian))) then
         ! Calculate the jacobian if requested (and either the numerical calculation is selected or the
         ! constitutive model does not provide the jacobian). Use approximation if strain increment is (almost) zero
         if (Norm(this%dot_strain) < setting_epsilon) then
            call this%Approximate_Jacobian(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_state=dot_state)
         else
            call this%Consistent_Jacobian(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_state=dot_state)
         end if
      end if
   end function Get_Dot_State
