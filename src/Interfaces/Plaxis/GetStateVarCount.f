! ------------------------------------------------------------------ !
subroutine GetStateVarCount(iModel, nStatevars)                      ! PLAXIS: Returns number of state variables of constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(out) :: nStatevars
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarCount     ! Export function name in dll
#endif

   nStatevars = 0
   if ((iModel >= 1) .and. (iModel <= size(numStateVariables))) then
      nStatevars = numStateVariables(iModel)
   end if
end subroutine GetStateVarCount
