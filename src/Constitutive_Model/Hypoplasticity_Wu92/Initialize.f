   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Hypoplasticity_Wu92), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)

      ! --- Parameters
      this%param_C1 = params(1)
      this%param_C2 = params(2)
      this%param_C3 = params(3)
      this%param_C4 = params(4)
   end subroutine Initialize
