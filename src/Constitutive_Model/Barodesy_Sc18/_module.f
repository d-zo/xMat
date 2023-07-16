! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Sc18_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Sc18
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_K_r, param_xi, param_kappa, param_e_c0, param_ce, &
                  param_c5, param_c6, param_c7
      real(dp) :: derived_c1, derived_c2, derived_c3, derived_c4, derived_c8, &
                  derived_zeta, derived_b_comp, derived_b_ext

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile function Calculate_Dot_State
end module Barodesy_Sc18_Class
