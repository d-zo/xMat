   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK23_Simpson), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n + 0.5h, y_n + 0.5h k_1)`
         cur_time=cur_time + 0.5_dp*dt, cur_state=cur_state + dt*0.5_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f(t_n + h, y_n - h k_1 + 2h k_2)`
         cur_time=cur_time + dt, cur_state=cur_state + dt*(-dot_state_k1 + 2.0_dp*dot_state_k2))

      new_state = cur_state + dt*(dot_state_k1/6.0_dp &              ! `y_{n+1} = y_n + \frac{1}{6}h k_1 + \frac{2}{3}h k_2 + \frac{1}{6}h k_3`
                + 2.0_dp/3.0_dp*dot_state_k2 + dot_state_k3/6.0_dp)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k3

      extra_state = cur_state + dt*dot_state_k2                      ! `y_{p,n+1} = y_n + h k_2`
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
