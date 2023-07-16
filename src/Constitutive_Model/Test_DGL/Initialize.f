   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculateJacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Test_DGL), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculateJacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculateJacobian=calculateJacobian, provideJacobian=.False.)

      ! --- Parameters
      this%param_select = int(params(1))                             ! Select, which test DGL should be used
   end subroutine Initialize
