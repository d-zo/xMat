   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Fehlberg), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.25h, y_n+0.25h k_1)`
         cur_time=cur_time + 0.25_dp*dt, &
         cur_state=cur_state + dt*0.25_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{8}h, y_n+\frac{3}{32}h k_1 + \frac{9}{32}h k_2\right)`
         cur_time=cur_time + 0.375_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/32.0_dp*dot_state_k1 + 9.0_dp/32.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{12}{13}h, y_n+\frac{1932}{2197}h k_1 - \frac{7200}{2197}h k_2 + \frac{7296}{2197}h k_3\right)`
         cur_time=cur_time + 12.0_dp/13.0_dp*dt, &
         cur_state=cur_state + dt*(1932.0_dp/2197.0_dp*dot_state_k1 &
            - 7200.0_dp/2197.0_dp*dot_state_k2 + 7296.0_dp/2197.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+h, y_n+\frac{439}{216}h k_1 - 8h k_2 + \frac{3680}{513}h k_3 - \frac{845}{4104}h k_4\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(439.0_dp/216.0_dp*dot_state_k1 - 8.0_dp*dot_state_k2 &
            + 3680.0_dp/513.0_dp*dot_state_k3 - 845.0_dp/4104.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+\frac{1}{2}h, y_n-\frac{8}{27}h k_1 + 2h k_2 - \frac{3544}{2565}h k_3 + \frac{1859}{4104}h k_4 - \frac{11}{40}h k_5\right)`
         cur_time=cur_time + 0.5_dp*dt, &
         cur_state=cur_state + dt*(-8.0_dp/27.0_dp*dot_state_k1 + 2.0_dp*dot_state_k2 &
            - 3544.0_dp/2565.0_dp*dot_state_k3 + 1859.0_dp/4104.0_dp*dot_state_k4 &
            - 11.0_dp/40.0_dp*dot_state_k5))

      new_state = cur_state + dt*(25.0_dp/216.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{25}{216}h k_1 + \frac{1408}{2565}h k_3 + \frac{2197}{4104}h k_4 - \frac{1}{5}h k_5`
                + 1408.0_dp/2565.0_dp*dot_state_k3 + 2197.0_dp/4104.0_dp*dot_state_k4 &
                - 0.2_dp*dot_state_k5)

      extra_state = cur_state + dt*(16.0_dp/135.0_dp*dot_state_k1 &  ! `y_{p,n+1} =  y_n + \frac{16}{135}h k_1 + \frac{6656}{12825}h k_3 + \frac{28561}{56430}h k_4 - \frac{9}{50}h k_5 + \frac{2}{55}h k_6`
                + 6656.0_dp/12825.0_dp*dot_state_k3 + 28561.0_dp/56430.0_dp*dot_state_k4 &
                - 9.0_dp/50.0_dp*dot_state_k5 + 2.0_dp/55.0_dp*dot_state_k6)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
