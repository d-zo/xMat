   ! --------------------------------------------------------------- !
   function Valid_State(this, cur_stress, cur_voidratio, e_d, e_c, e_i, dt, dot_voidratio)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan, Write_Warning, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Nonzero_Division, Trace, Value_In_Interval
      !
      class(Hypoplasticity_VW96), intent(in) :: this
      real(dp), dimension(__matrix__), intent(in) :: cur_stress
      real(dp), intent(inout) :: cur_voidratio
      real(dp), intent(out) :: e_d, e_c, e_i
      real(dp), intent(in) :: dt
      real(dp), intent(out) :: dot_voidratio
      logical :: Valid_State
      ! ------------------------------------------------------------ !
      real(dp) :: trT, Bauer_Term, temp_voidratio

      Valid_State = .False.
      dot_voidratio = 0.0_dp
      temp_voidratio = cur_voidratio

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      if (Is_Nan(trT)) then
         call Write_Error_And_Exit('Valid_State: trT is NaN')
      end if

      if (this%Is_Valid_Stress_State(cur_stress)) then
         ! Compute the current void ratios according to Bauer (1996)
         Bauer_Term = exp(-(-trT/this%param_h_s)**this%param_n_H)    ! Exponential compression law by Bauer (1996): `B(\tr{(\mathbf{T})}) = \exp{}[-(\frac{\tr{(\mathbf{T})}}{h_s})^n_H]`
         e_d = this%param_e_d0 * Bauer_Term                          ! Current densest void ratio  `e_d = e_{d0} B(\tr{(\mathbf{T})})`
         e_c = this%param_e_c0 * Bauer_Term                          ! Current critical void ratio `e_c = e_{c0} B(\tr{(\mathbf{T})})`
         e_i = this%param_e_i0 * Bauer_Term                          ! Current loosest void ratio  `e_i = e_{i0} B(\tr{(\mathbf{T})})`

         ! Assuming void ratio is not specified as statevariable try materialparameter
         if (cur_voidratio < setting_epsilon) then
            cur_voidratio = this%param_e_ini * Bauer_Term
         end if

         if (.not. Value_In_Interval(cur_voidratio, [e_d, e_i])) then
            ! Only give feedback if it is considerably off
            if (.not. Value_In_Interval(cur_voidratio, [0.96_dp*e_d, 1.04_dp*e_i])) then
               call Write_Warning('Valid_State: ' // Formatval('void ratio:', cur_voidratio) // &
                  ' out of interval [' // Formatval('e_d:', e_d) // ', ' // Formatval('e_i:', e_i) // ']')
            end if

            cur_voidratio = max(e_d, min(e_i, cur_voidratio))
         end if
         Valid_State = .True.
      end if

      dot_voidratio = Nonzero_Division(val=(cur_voidratio - temp_voidratio), fac=dt)
   end function Valid_State
