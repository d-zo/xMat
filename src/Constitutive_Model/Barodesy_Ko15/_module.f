! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Ko15_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Ko15
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_c2, param_c3, param_c4, param_c5, param_e_c0, param_e_min
      real(dp) :: derived_c1

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile function Calculate_Dot_State
end module Barodesy_Ko15_Class
