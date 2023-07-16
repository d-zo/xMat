   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Dormand_Prince), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, dot_state_k7, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.2h, y_n+0.2h k_1)`
         cur_time=cur_time + 0.2_dp*dt, &
         cur_state=cur_state + dt*0.2_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{10}h, y_n+\frac{3}{40}h k_1 + \frac{9}{40}h k_2\right)`
         cur_time=cur_time + 0.3_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/40.0_dp*dot_state_k1 + 9.0_dp/40.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{4}{5}h, y_n+\frac{44}{45}h k_1 - \frac{56}{15}h k_2 + \frac{32}{9}h k_3\right)`
         cur_time=cur_time + 0.8_dp*dt, &
         cur_state=cur_state + dt*(44.0_dp/45.0_dp*dot_state_k1 - 56.0_dp/15.0_dp*dot_state_k2 &
            + 32.0_dp/9.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+\frac{8}{9}h, y_n+\frac{19372}{6561}h k_1 - \frac{25360}{2187}h k_2 + \frac{64448}{6561}h k_3 - \frac{212}{729}h k_4\right)`
         cur_time=cur_time + 8.0_dp/9.0_dp*dt, &
         cur_state=cur_state + dt*(19372.0_dp/6561.0_dp*dot_state_k1 &
            - 25360.0_dp/2187.0_dp*dot_state_k2 + 64448.0_dp/6561.0_dp*dot_state_k3 &
            - 212.0_dp/729.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+h, y_n+\frac{9017}{3168}h k_1 - \frac{355}{33}h k_2 + \frac{46732}{5247}h k_3 + \frac{49}{176}h k_4 - \frac{5103}{18656}h k_5\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(9017.0_dp/3168.0_dp*dot_state_k1 - 355.0_dp/33.0_dp*dot_state_k2 &
            + 46732.0_dp/5247.0_dp*dot_state_k3 + 49.0_dp/176.0_dp*dot_state_k4 &
            - 5103.0_dp/18656.0_dp*dot_state_k5))

      new_state = cur_state + dt*(35.0_dp/384.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{35}{384}h k_1 + \frac{500}{1113}h k_3 + \frac{125}{192}h k_4 - \frac{2187}{6784}h k_5 + \frac{11}{84}h k_6`
                + 500.0_dp/1113.0_dp*dot_state_k3 + 125.0_dp/192.0_dp*dot_state_k4 &
                - 2187.0_dp/6784.0_dp*dot_state_k5 + 11.0_dp/84.0_dp*dot_state_k6)

      dot_state_k7 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_7 = f(t_n+h, y_{p,n+1})`
         cur_time=cur_time + dt, &
         cur_state=new_state)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k7

      extra_state = cur_state &                                      ! `y_{p,n+1} =  y_n + \frac{5179}{57600}h k_1 + \frac{7571}{16695}h k_3 + \frac{393}{640}h k_4 - \frac{92097}{339200}h k_5 + \frac{187}{2100}h k_6 + \frac{1}{40}h k_7`
                  + dt*(5179.0_dp/57600.0_dp*dot_state_k1 + 7571.0_dp/16695.0_dp*dot_state_k3 &
                  + 393.0_dp/640.0_dp*dot_state_k4 - 92097.0_dp/339200.0_dp*dot_state_k5 &
                  + 187.0_dp/2100.0_dp*dot_state_k6 + 1.0_dp/40.0_dp*dot_state_k7)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
