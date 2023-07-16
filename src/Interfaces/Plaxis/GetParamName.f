! ------------------------------------------------------------------ !
subroutine GetParamName(iModel, iParameter, ParameterName)           ! PLAXIS: Returns name of selected parameter of constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelParameterNamesUnits, numParameters
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
   integer, intent(in) :: iParameter
   character(len=*), intent(out) :: ParameterName
   ! --------------------------------------------------------------- !
   integer :: lenName, num_params
   character(len=16) :: tempParameterName

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetParamName         ! Export function name in dll
#endif

   tempParameterName = '                '
   if ((iModel >= 1) .and. (iModel <= size(ModelParameterNamesUnits, 2))) then
      ! List not only parameters but also state variables to be able to initialize the latter by additional values in the former
      num_params = numParameters(iModel)
      if (iParameter > num_params) then
         call GetStateVarNameInternal(iModel=iModel, iStatevar=iParameter-num_params, StatevarName=tempParameterName)
         tempParameterName = tempParameterName(2:10)
         if (tempParameterName /= '') then
            tempParameterName = 'state: ' // tempParameterName
         end if
      else
         if ((iParameter >= 1) .and. (iParameter <= numParameters(iModel)) &
            .and. (iParameter <= size(ModelParameterNamesUnits, 1))) then

            ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
            tempParameterName = ModelParameterNamesUnits(iParameter, iModel)
            tempParameterName(10:16) = '       '
         end if
      end if
   end if

   lenName = len_trim(tempParameterName)
   ! The first char represents the binary length of the string
   ParameterName = char(lenName) // tempParameterName(1:lenName)
end subroutine GetParamName
