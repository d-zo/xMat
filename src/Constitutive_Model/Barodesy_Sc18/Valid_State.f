   ! --------------------------------------------------------------- !
   function Valid_State(this, cur_stress)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Is_Nan, Write_Error_And_Exit
      use Math_Operations, only: Trace
      !
      class(Barodesy_Sc18), intent(in) :: this
      real(dp), dimension(__matrix__), intent(in) :: cur_stress
      logical :: Valid_State
      ! ------------------------------------------------------------ !
      real(dp) :: trT

      Valid_State = .False.

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      if (Is_Nan(trT)) then
         call Write_Error_And_Exit('Valid_State: trT is NaN')
      end if

      if (this%Is_Valid_Stress_State(cur_stress)) then
         Valid_State = .True.
      end if
   end function Valid_State
