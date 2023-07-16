! ==================================================================================================================== !
module Solver_Baseclass
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model

   implicit none

   private
   public :: Solver

   ! --------------------------------------------------------------- !
   type, abstract :: Solver
   ! --------------------------------------------------------------- !
      character(len=24) :: name
      logical :: current_error_assigned, next_dot_state_assigned, stepsize_fixed, last_integ_success
      integer :: last_num_steps_accepted, last_num_steps_rejected
      real(dp), dimension(setting_max_internal_states) :: current_error, next_dot_state
      logical, dimension(setting_max_internal_states) :: mask_ignored_indices
      real(dp) :: exp_stepgrow, max_stepgrow, exp_stepshrink, min_stepshrink

      contains

      procedure(integration_algorithm_interface), deferred :: Integration_Algorithm

      procedure :: Initialize
      procedure :: Integrate
      procedure :: Set_State_Mask
      procedure :: Set_Error_Estimate
      procedure :: Get_Statistics
   end type

   ! --------------------------------------------------------------- !
   abstract interface
   ! --------------------------------------------------------------- !
      function integration_algorithm_interface(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) &
         result(new_state)
         import
         class(Solver), intent(inout) :: this
         class(Constitutive_Model), intent(inout) :: integ_obj
         real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
         real(dp), intent(in) :: cur_time, dt
         real(dp), dimension(setting_max_internal_states) :: new_state
      end function
   end interface


   contains


#addfile subroutine Initialize


#addfile subroutine Integrate


#addfile function Set_State_Mask


#addfile function Set_Error_Estimate


#addfile function Get_Statistics
end module Solver_Baseclass
