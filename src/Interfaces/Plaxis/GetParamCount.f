! ------------------------------------------------------------------ !
subroutine GetParamCount(iModel, nParameters)                        ! PLAXIS: Returns number of parameters of chosen constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: numParameters, numStatevariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(out) :: nParameters
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetParamCount        ! Export function name in dll
#endif

   nParameters = 0
   if ((iModel >= 1) .and. (iModel <= size(numParameters))) then
      ! Return number of parameters and statevariables to be able to initialize the latter by additional values in the former
      nParameters = numParameters(iModel) + numStatevariables(iModel)
   end if
end subroutine GetParamCount
