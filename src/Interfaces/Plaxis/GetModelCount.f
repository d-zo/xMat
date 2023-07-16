! ------------------------------------------------------------------ !
subroutine GetModelCount(numModels)                                  ! PLAXIS: Returns amount of usable constitutive models
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: numParameters
   implicit none
   !
   integer, intent(out) :: numModels
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetModelCount        ! Export function name in dll
#endif

   numModels = size(numParameters)
end subroutine GetModelCount
