! ==================================================================================================================== !
module Xmat_Class
   use General_Settings, only: dp, setting_len_id, setting_num_statevariables, setting_max_internal_states
   use Solver_Class, only: Solver, Select_Solver
   use Constitutive_Model_Class, only: Constitutive_Model, Select_Constitutive_Model

   implicit none

   private
   public :: Xmat, Xmat_Initialize

   ! --------------------------------------------------------------- !
   type :: Xmat
   ! --------------------------------------------------------------- !
      private
      class(Solver), allocatable :: xmat_solver
      class(Constitutive_Model), allocatable :: xmat_constitutive_model

      real(dp), dimension(setting_num_statevariables) :: statevariables
      real(dp) :: assigned_time, assigned_dt
      real(dp), dimension(__matrix__) ::  assigned_dot_strain, saved_stress
      integer :: nstates

      contains

      procedure :: Import_State
      procedure :: Calculate_Results
      procedure :: Dot_Results
   end type


   contains


#addfile subroutine Xmat_Initialize


#addfile subroutine Import_State


#addfile subroutine Calculate_Results


#addfile subroutine Dot_Results
end module Xmat_Class
