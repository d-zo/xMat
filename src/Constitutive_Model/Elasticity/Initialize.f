   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval
      !
      class(Elasticity), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)

      ! --- Parameters and description
      this%param_youngs_modulus = params(1)                          ! Young's modulus `E` (in kPa)
      this%param_nu             = params(2)                          ! Poisson's ratio `\nu`

      if (firstcall) then                                            ! Upper limit on Young's modulus is chosen arbitrarily
         call Abort_If_Not_In_Interval('E', this%param_youngs_modulus, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
      end if
   end subroutine Initialize
