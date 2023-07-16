! ------------------------------------------------------------------ !
pure subroutine GetStateVarUnitInternal(iModel, iStatevar, StatevarUnit)
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarUnit
   ! --------------------------------------------------------------- !
   integer :: lenUnit
   character(len=16) :: tempStatevarUnit

   tempStatevarUnit = '-               '
   if ((iModel >= 1) .and. (iModel <= size(ModelStatevarNamesUnits, 2))) then
      if ((iStatevar >= 1) .and. (iStatevar <= numStateVariables(iModel)) &
         .and. (iStatevar <= size(ModelStatevarNamesUnits, 1))) then

         ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
         tempStatevarUnit = ModelStatevarNamesUnits(iStatevar, iModel)
         tempStatevarUnit = tempStatevarUnit(11:16) // '          '
      end if
   end if
   lenUnit = len_trim(tempStatevarUnit)
   ! The first char represents the binary length of the string
   StatevarUnit = char(lenUnit) // tempStatevarUnit(1:lenUnit)
end subroutine GetStateVarUnitInternal
