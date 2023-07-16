   ! --------------------------------------------------------------- !
   function Valid_State(this, cur_stress, cur_voidratio, pressure_equiv)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan, Write_Error_And_Exit
      use Math_Operations, only: Trace
      !
      class(Viscohypoplasticity_Ni03), intent(in) :: this
      real(dp), dimension(__matrix__), intent(in) :: cur_stress
      real(dp), intent(inout) :: cur_voidratio
      real(dp), intent(out) :: pressure_equiv
      logical :: Valid_State
      ! ------------------------------------------------------------ !
      real(dp) :: trT

      Valid_State = .False.

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      if (Is_Nan(trT)) then
         call Write_Error_And_Exit('Valid_State: trT is NaN')
      end if

      if (cur_voidratio < setting_epsilon) then
         call Write_Error_And_Exit('Valid_State: void ratio invalid (too small) and no default value available')
      end if

      if (this%Is_Valid_Stress_State(cur_stress)) then
         pressure_equiv = 100.0_dp*((1.0_dp + this%param_e_100) &    ! Calculate transformed (4.44) of Niemunis (2003): `T_e = T_{e0}\left(\frac{1+e}{1+e_{100}}\right)^{-\frac{1}{\lambda}}`
                        / (1.0_dp + cur_voidratio)) &                ! with `T_{e0} = 100`
                        **(1.0_dp/this%param_lambda)
         Valid_State = .True.
      end if
   end function Valid_State
