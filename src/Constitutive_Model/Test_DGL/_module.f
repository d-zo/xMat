! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Test_DGL_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Test_DGL
   ! --------------------------------------------------------------- !
      private
      integer :: param_select

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


#addfile subroutine Initialize


#addfile function Calculate_Dot_State
end module Test_DGL_Class
