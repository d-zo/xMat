! ------------------------------------------------------------------ !
subroutine GetStateVarUnit(iModel, iStatevar, StatevarUnit)          ! PLAXIS: Returns unit of state variable of constitutive model
! ------------------------------------------------------------------ !
   implicit none
   !
   interface
      subroutine GetStateVarUnitInternal(iModel, iStatevar, StatevarUnit)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarUnit
      end subroutine GetStateVarUnitInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarUnit
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarUnit      ! Export function name in dll
#endif

   call GetStateVarUnitInternal(iModel=iModel, iStatevar=iStatevar, StatevarUnit=StatevarUnit)
end subroutine GetStateVarUnit
