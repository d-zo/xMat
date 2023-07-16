   ! --------------------------------------------------------------- !
   function Get_Pressure_Preconsol(this, overcrit, trT)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      class(Viscohypoplasticity_Ni03), intent(in) :: this
      real(dp), intent(in) :: overcrit, trT
      real(dp) :: Get_Pressure_Preconsol
      ! ------------------------------------------------------------ !
      real(dp) :: p_mean

      p_mean = -trT/3.0_dp
      if (overcrit > 1.0_dp) then
         Get_Pressure_Preconsol = p_mean/2.0_dp*(1.0_dp + overcrit)*(1.0_dp + this%param_beta_b)
      else
         Get_Pressure_Preconsol = p_mean*(1.0_dp - this%param_beta_b &
                                * sqrt(1.0_dp + overcrit*(this%param_beta_b**2 - 1.0_dp)))/(1.0_dp - this%param_beta_b)
      end if

      if (Get_Pressure_Preconsol < 0.0_dp) then
         call Write_Error_And_Exit('Get_Pressure_Preconsol: Negative preconsolidation pressure ' // &
            Formatval('p_e^+: ', Get_Pressure_Preconsol))
      end if
   end function Get_Pressure_Preconsol
