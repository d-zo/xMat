   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK23_Bogacki_Shampine), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.5h, y_n+0.5h k_1)`
         cur_time=cur_time + 0.5_dp*dt, cur_state=cur_state + dt*0.5_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f(t_n+0.75h, y_n+0.75h k_2)`
         cur_time=cur_time + 0.75_dp*dt, cur_state=cur_state + dt*0.75_dp*dot_state_k2)

      new_state = cur_state + dt*(2.0_dp/9.0_dp*dot_state_k1 &       ! `y_{n+1} = y_n + \frac{2}{9}h k_1 + \frac{1}{3}h k_2 + \frac{4}{9}h k_3`
                + 1.0_dp/3.0_dp*dot_state_k2 + 4.0_dp/9.0_dp*dot_state_k3)

      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f(t_n+h, y_{n+1})`
         cur_time=cur_time + dt, cur_state=new_state)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k4

      extra_state = cur_state + dt*(7.0_dp/24.0_dp*dot_state_k1 &    ! `y_{p,n+1} = y_n + \frac{7}{24}h k_1 + \frac{1}{4}h k_2 + \frac{1}{3}h k_3 + \frac{1}{8}h k_4`
                  + 1.0_dp/4.0_dp*dot_state_k2 + 1.0_dp/3.0_dp*dot_state_k3 + 1.0_dp/8.0_dp*dot_state_k4)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
