! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Hypoplasticity_VW96_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Hypoplasticity_VW96
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_nu, param_h_s, param_n_H, param_e_d0, param_e_c0, param_e_i0, &
                  param_alpha_H, param_beta_H, param_m_T, param_m_R, param_R_max, param_beta_R, param_chi, &
                  param_e_ini
      real(dp) :: derived_a, derived_a2, derived_f_b_lastterm, derived_b2

      contains

      procedure :: Initialize
      procedure :: Valid_State
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile subroutine Valid_State


#addfile function Calculate_Dot_State
end module Hypoplasticity_VW96_Class
