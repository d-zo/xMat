! ------------------------------------------------------------------ !
subroutine GetParamUnit(iModel, iParameter, ParameterUnit)           ! PLAXIS: Returns unit of parameter of chosen constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelParameterNamesUnits, numParameters
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
   integer, intent(in) :: iParameter
   character(len=*) :: ParameterUnit
   ! --------------------------------------------------------------- !
   integer :: lenUnit, num_params
   character(len=16) :: tempParameterUnit

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetParamUnit         ! Export function name in dll
#endif

   tempParameterUnit = '-               '
   if ((iModel >= 1) .and. (iModel <= size(ModelParameterNamesUnits, 2))) then
      ! List not only parameters but also state variables to be able to initialize the latter by additional values in the former
      num_params = numParameters(iModel)
      if (iParameter > num_params) then
         call GetStateVarUnitInternal(iModel=iModel, iStatevar=iParameter-num_params, StatevarUnit=tempParameterUnit)
         tempParameterUnit = tempParameterUnit(2:16) // ' '
      else
         if ((iParameter >= 1) .and. (iParameter <= numParameters(iModel)) &
            .and. (iParameter <= size(ModelParameterNamesUnits, 1))) then

            ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
            tempParameterUnit = ModelParameterNamesUnits(iParameter, iModel)
            tempParameterUnit = tempParameterUnit(11:16) // '          '
         end if
      end if
   end if
   lenUnit = len_trim(tempParameterUnit)
   ! The first char represents the binary length of the string
   ParameterUnit = char(lenUnit) // tempParameterUnit(1:lenUnit)
end subroutine GetParamUnit
