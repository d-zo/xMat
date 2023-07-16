! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK45_Cash_Karp_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK45_Cash_Karp
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


#addfile function Integration_Algorithm
end module RK45_Cash_Karp_Class
