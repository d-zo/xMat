! ------------------------------------------------------------------ !
subroutine GetStateVarName(iModel, iStatevar, StatevarName)          ! PLAXIS: Returns name of state variable of chosen constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   interface
      subroutine GetStateVarNameInternal(iModel, iStatevar, StatevarName)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarName
      end subroutine GetStateVarNameInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarName
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarName      ! Export function name in dll
#endif

   call GetStateVarNameInternal(iModel=iModel, iStatevar=iStatevar, StatevarName=StatevarName)
end subroutine GetStateVarName
