! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Elasticity_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Elasticity
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_youngs_modulus, param_nu

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile function Calculate_Dot_State
end module Elasticity_Class
