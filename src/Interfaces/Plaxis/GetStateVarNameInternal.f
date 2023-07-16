! ------------------------------------------------------------------ !
subroutine GetStateVarNameInternal(iModel, iStatevar, StatevarName)
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarName
   ! --------------------------------------------------------------- !
   integer :: lenName
   character(len=16) :: tempStatevarName

   tempStatevarName = '                '
   if ((iModel >= 1) .and. (iModel <= size(ModelStatevarNamesUnits, 2))) then
      if ((iStatevar >= 1) .and. (iStatevar <= numStateVariables(iModel)) &
         .and. (iStatevar <= size(ModelStatevarNamesUnits, 1))) then

         ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
         tempStatevarName = ModelStatevarNamesUnits(iStatevar, iModel)
         tempStatevarName(10:16) = '       '
      end if
   end if
   lenName = len_trim(tempStatevarName)
   ! The first char represents the binary length of the string
   StatevarName = char(lenName) // tempStatevarName(1:lenName)
end subroutine GetStateVarNameInternal
