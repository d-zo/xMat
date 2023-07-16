! ==================================================================================================================== !
module Constitutive_Model_Baseclass
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   implicit none

   private
   public :: Constitutive_Model

   ! --------------------------------------------------------------- !
   type, abstract :: Constitutive_Model
   ! --------------------------------------------------------------- !
      real(dp) :: overall_dt
      real(dp), dimension(__matrix__) :: dot_strain
      logical :: calculateJacobian
      logical :: provideJacobian
      real(dp), dimension(setting_num_statevariables) :: direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask

      contains

      procedure(initialize_interface), deferred :: Initialize
      procedure(calculate_dot_state_interface), deferred :: Calculate_Dot_State

      procedure :: Base_Initialization
      procedure :: Get_Dot_State
      procedure :: Set_Values
      procedure :: Get_Direct_Variables
      procedure :: Get_Direct_Variables_Mask
      procedure :: Approximate_Jacobian
      procedure :: Consistent_Jacobian
      procedure :: Elasticity
      procedure, nopass :: Is_Valid_Stress_State
      procedure, nopass :: Get_Factor_F
      procedure :: Intergranular_Strain
   end type

   ! --------------------------------------------------------------- !
   abstract interface
   ! --------------------------------------------------------------- !
      subroutine initialize_interface(this, params, calculateJacobian, firstcall)
         import
         class(Constitutive_Model), intent(inout) :: this
         real(dp), dimension(:), intent(in) :: params
         logical, intent(in) :: calculateJacobian, firstcall
      end subroutine initialize_interface

      function calculate_dot_state_interface(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
         import
         class(Constitutive_Model), intent(inout) :: this
         real(dp), intent(in) :: ref_dt
         real(dp), intent(in) :: cur_time
         real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
         real(dp), dimension(__matrix__), intent(in) :: dot_strain
         real(dp), dimension(setting_max_internal_states) :: dot_state
      end function calculate_dot_state_interface
   end interface


   contains


#addfile subroutine Base_Initialization


#addfile function Get_Dot_State


#addfile subroutine Set_Values


#addfile subroutine Get_Direct_Variables


#addfile function Get_Direct_Variables_Mask


#addfile subroutine Approximate_Jacobian


#addfile subroutine Consistent_Jacobian


#addfile subroutine Elasticity


#addfile function Is_Valid_Stress_State


#addfile function Get_Factor_F


#addfile subroutine Intergranular_Strain
end module Constitutive_Model_Baseclass
