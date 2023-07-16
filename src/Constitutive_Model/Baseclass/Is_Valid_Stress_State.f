   ! --------------------------------------------------------------- !
   pure function Is_Valid_Stress_State(stress) result(in_limits)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_very_small_stress
      use Math_Operations, only: Trace
      !
      real(dp), dimension(__matrix__), intent(in) :: stress
      logical :: in_limits
      ! ------------------------------------------------------------ !
      real(dp) :: trT

      ! NOTE: Currently this function checks the trace only. For a more reliable check the individual diagonal components should be checked
      !       and possibly also the actually used limit criterion (Mohr-Coulomb, Matsuoka-Nakai, ...)
      in_limits = .True.
      trT = Trace(stress)
      if (trT > setting_very_small_stress) then
         in_limits = .False.
      end if
   end function Is_Valid_Stress_State
