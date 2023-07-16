   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Cash_Karp), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.2h, y_n+0.2h k_1)`
         cur_time=cur_time + 0.2_dp*dt, &
         cur_state=cur_state + dt*0.2_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{10}h, y_n+\frac{3}{40}h k_1 + \frac{9}{40}h k_2\right)`
         cur_time=cur_time + 0.3_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/40.0_dp*dot_state_k1 + 9.0_dp/40.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{3}{5}h, y_n+\frac{3}{10}h k_1 - \frac{9}{10}h k_2 + \frac{6}{5}h k_3\right)`
         cur_time=cur_time + 0.6_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/10.0_dp*dot_state_k1 - 9.0_dp/10.0_dp*dot_state_k2 &
            + 6.0_dp/5.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+h, y_n-\frac{11}{54}h k_1 + \frac{5}{2}h k_2 - \frac{70}{27}h k_3 + \frac{35}{27}h k_4\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(-11.0_dp/54.0_dp*dot_state_k1 + 5.0_dp/2.0_dp*dot_state_k2 &
            - 70.0_dp/27.0_dp*dot_state_k3 + 35.0_dp/27.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+\frac{7}{8}h, y_n+\frac{1631}{55296}h k_1 + \frac{175}{512}h k_2 + \frac{575}{13824}h k_3 + \frac{44275}{110592}h k_4 + \frac{253}{4096}h k_5\right)`
         cur_time=cur_time + 7.0_dp/8.0_dp*dt, &
         cur_state=cur_state + dt*(1631.0_dp/55296.0_dp*dot_state_k1 &
            + 175.0_dp/512.0_dp*dot_state_k2 + 575.0_dp/13824.0_dp*dot_state_k3 &
            + 44275.0_dp/110592.0_dp*dot_state_k4 + 253.0_dp/4096.0_dp*dot_state_k5))

      new_state = cur_state + dt*(37.0_dp/378.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{37}{378}h k_1 + \frac{250}{621}h k_3 + \frac{125}{594}h k_4 - \frac{512}{1771}h k_6`
                + 250.0_dp/621.0_dp*dot_state_k3 + 125.0_dp/594.0_dp*dot_state_k4 &
                + 512.0_dp/1771.0_dp*dot_state_k6)

      extra_state = cur_state &                                      ! `y_{p,n+1} =  y_n + \frac{2825}{27648}h k_1 + \frac{18575}{48384}h k_3 + \frac{13525}{55296}h k_4 + \frac{277}{14336}h k_5 + \frac{1}{4}h k_6`
                  + dt*(2825.0_dp/27648.0_dp*dot_state_k1 + 18575.0_dp/48384.0_dp*dot_state_k3 &
                  + 13525.0_dp/55296.0_dp*dot_state_k4 + 277.0_dp/14336.0_dp*dot_state_k5 + 0.25_dp*dot_state_k6)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
