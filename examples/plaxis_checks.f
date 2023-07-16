program plaxis_check
   implicit none

   integer :: numModels, iModel, iParameter, nParameters, iStatevar, nStatevars
   character(len=30) :: ModelName, ParameterName, ParameterUnit, StatevarName, StatevarUnit

   write(*, '(a,a,a)') char(10), '# Checking GetModelCount', char(10)
   call GetModelCount(numModels)
   write(*, '(i1)') numModels

   write(*, '(a,a,a)') char(10), '# Checking GetModelName', char(10)
   do iModel = 0, 9
      ModelName = repeat('.', 30)
      call GetModelName(iModel, ModelName)
      write(*, '(i1,a1,a)') iModel, ' ', ModelName(2:30)
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetParamCount', char(10)
   do iModel = 0, 9
      nParameters = 0
      call GetParamCount(iModel, nParameters)
      write(*, '(i1,a1,i2)') iModel, ' ', nParameters
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetParamName', char(10)
   do iModel = 0, 7
      do iParameter = 0, 17
         ParameterName = repeat('.', 30)
         call GetParamName(iModel, iParameter, ParameterName)
         write(*, '(i1,a2,i2,a1,a)') iModel, ', ', iParameter, ' ', ParameterName(2:30)
      end do
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetParamUnit', char(10)
   do iModel = 0, 7
      do iParameter = 0, 17
         ParameterUnit = repeat('.', 30)
         call GetParamUnit(iModel, iParameter, ParameterUnit)
         write(*, '(i1,a2,i2,a1,a)') iModel, ', ', iParameter, ' ', ParameterUnit(2:30)
      end do
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetStatevarCount', char(10)
   do iModel = 0, 9
      nStatevars = 0
      call GetStatevarCount(iModel, nStatevars)
      write(*, '(i1,a1,i2)') iModel, ' ', nStatevars
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetStatevarName', char(10)
   do iModel = 0, 7
      do iStatevar = 0, 12
         StatevarName = repeat('.', 30)
         call GetStatevarName(iModel, iStatevar, StatevarName)
         write(*, '(i1,a2,i2,a1,a)') iModel, ', ', iStatevar, ' ', StatevarName(2:30)
      end do
   end do

   write(*, '(a,a,a)') char(10), '# Checking GetStatevarUnit', char(10)
   do iModel = 0, 7
      do iStatevar = 0, 12
         StatevarUnit = repeat('.', 30)
         call GetStatevarUnit(iModel, iStatevar, StatevarUnit)
         write(*, '(i1,a2,i2,a1,a)') iModel, ', ', iStatevar, ' ', StatevarUnit(2:30)
      end do
   end do
end program plaxis_check
