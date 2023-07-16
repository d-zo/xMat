   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculateJacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Hypoplasticity_Wu92), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculateJacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculateJacobian=calculateJacobian, provideJacobian=.True.)

      ! --- Parameters
      this%param_C1 = params(1)
      this%param_C2 = params(2)
      this%param_C3 = params(3)
      this%param_C4 = params(4)
   end subroutine Initialize
