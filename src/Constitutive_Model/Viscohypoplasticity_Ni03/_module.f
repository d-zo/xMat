! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Viscohypoplasticity_Ni03_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Viscohypoplasticity_Ni03
   ! --------------------------------------------------------------- !
      private
      ! NOTE: Currently it is confusing to have four very similar named but different parameters:
      !         - param_beta_b (viscohypoplastic parameter)
      !         - param_beta_R (extended hypoplastic parameter)
      !         - derived_beta_b (internal parameter)
      !         - derived_b2 (internal parameter)

      real(dp) :: param_e_100, param_nu, param_lambda, param_kappa, param_beta_b, param_I_v, param_D_r, &
                  param_phi_c, param_m_T, param_m_R, param_R_max, param_beta_R, param_Chi, param_OCR
      real(dp) :: derived_a, derived_a2, derived_pq_incl, derived_beta_b, derived_b2

      contains

      procedure :: Initialize
      procedure :: Valid_State
      procedure :: Get_Pressure_Preconsol
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile subroutine Valid_State


#addfile subroutine Get_Pressure_Preconsol


#addfile function Calculate_Dot_State
end module Viscohypoplasticity_Ni03_Class
