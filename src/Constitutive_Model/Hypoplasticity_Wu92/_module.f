! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Hypoplasticity_Wu92_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Hypoplasticity_Wu92
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_C1, param_C2, param_C3, param_C4

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile function Calculate_Dot_State
end module Hypoplasticity_Wu92_Class
