! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Ko21_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Ko21
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_c2, param_c3, param_c4, param_c5, param_kappa1, param_kappa2, &
                  param_k1, param_k2, param_e_c0
      real(dp) :: derived_c1

      contains

      procedure :: Initialize
      procedure :: Valid_State
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile subroutine Valid_State


#addfile function Calculate_Dot_State
end module Barodesy_Ko21_Class
