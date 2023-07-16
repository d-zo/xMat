! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK23_Simpson_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK23_Simpson
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


#addfile function Integration_Algorithm
end module RK23_Simpson_Class
