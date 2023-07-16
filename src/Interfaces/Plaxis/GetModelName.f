! ------------------------------------------------------------------ !
subroutine GetModelName(iModel, ModelName)                           ! PLAXIS: Returns name of usable constitutive models
! ------------------------------------------------------------------ !
   use General_Settings, only: setting_len_id
   use PlaxisInformationPool, only: ModelList
   implicit none
   !
   integer, intent(in) :: iModel
   character(len=*), intent(out) :: ModelName
   ! --------------------------------------------------------------- !
   character(len=setting_len_id) :: tempModelName

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetModelName         ! Export function name in dll
#endif

   tempModelName = repeat('-', setting_len_id)
   if ((iModel >= 1) .and. (iModel <= size(ModelList))) then
      tempModelName = ModelList(iModel)
   end if
   ! The first char represents the binary length of the string
   ModelName = char(setting_len_id) // tempModelName(1:setting_len_id)
end subroutine GetModelName
