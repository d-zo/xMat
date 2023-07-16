   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(Euler_Explicit), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      new_state = cur_state + dt*dot_state_k1                        ! `y_{n+1} = y_n + h f(t_n, y_n)`
   end function Integration_Algorithm
