! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                       m e x  F u n c t i o n                       ! Entry point for Matlab
! ------------------------------------------------------------------ !
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
! ------------------------------------------------------------------ !
   use General_Settings, only: dp, setting_len_id, setting_solver_default, setting_max_mat_params, &
                               setting_num_statevariables, Uppercase, Estimate_Components, Check_Input_Dimensions, &
                               Write_Error_And_Exit
   use Debug, only: Formatval
   use Math_Operations, only: Nonzero_Division, Import_Matrix, Export_Matrix, Export_Tensor
   use Xmat_Class, only: Xmat, Xmat_Initialize
   implicit none
   !
   integer(kind=4) :: nlhs, nrhs
   mwpointer :: plhs(nlhs), prhs(nrhs)
   ! --------------------------------------------------------------- !
   ! --- Auxiliary input/output variables
   mwpointer :: mxGetNumberOfElements, mxGetData, mxCreateDoubleMatrix, mxGetPr
   real(kind=8) :: dt, tt
   real(kind=8), dimension(6) :: stress, strain                      ! Assume Matlab input to be vec6
   real(kind=8), dimension(6, 6) :: jacobian
   real(dp), dimension(6) :: temp_mat, tmp_dot_strain
   mwSize :: nel_matname, nel_props, nel_statev, nel_stress, nel_strain, nel_dt, nel_tt
   real(kind=8), dimension(setting_max_mat_params) :: matprops
   real(kind=8), dimension(setting_num_statevariables) :: state
   ! --- Variables with custom precision for the actual calculations
   integer :: num_dimensions, num_shear
   real(dp), dimension(setting_max_mat_params) :: materialproperties
   real(dp), dimension(setting_num_statevariables) :: inp_state
   real(dp), dimension(__matrix__) :: inp_stress, inp_strain, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt
   character(len=setting_len_id) :: matname

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(__tensor__) :: jacobian_internal
   character(len=setting_len_id) :: identifier
   character(len=13), parameter :: intername = 'mex_function '
   integer :: idx
   logical :: firstcall

   ! The function should be called in matlab with the following seven arguments (in order):
   ! materialname, materialproperties, inp_state, inp_stress, inp_dot_strain, timeincrement, totaltime
   if (nrhs .ne. 7) then
      call Write_Error_And_Exit('mexFunction: Expecting exactly 7 input arguments from Matlab')
   end if
   ! The function returns the following two arguments to matlab (in order):
   ! newstress, newstate, jacobian
   if (nlhs .ne. 3) then
     call Write_Error_And_Exit('mexFunction: Providing exactly 3 outputs to Matlab')
   end if

   ! --- Import variables from Matlab pointers
   nel_matname = mxGetNumberOfElements(prhs(1))                      ! If nel_matname > len(matname), all extra characters are ignored
   nel_props   = mxGetNumberOfElements(prhs(2))
   nel_statev  = mxGetNumberOfElements(prhs(3))
   nel_stress  = mxGetNumberOfElements(prhs(4))
   nel_strain  = mxGetNumberOfElements(prhs(5))
   nel_dt      = mxGetNumberOfElements(prhs(6))
   nel_tt      = mxGetNumberOfElements(prhs(7))

   ! NOTE: There is currently no check that the number of elements match for each vector
   call mxGetString(prhs(1), matname, nel_matname)
   call mxCopyPtrToReal8(mxGetData(prhs(2)), matprops, nel_props)
   call mxCopyPtrToReal8(mxGetData(prhs(3)), state, nel_statev)
   call mxCopyPtrToReal8(mxGetData(prhs(4)), stress, nel_stress)
   call mxCopyPtrToReal8(mxGetData(prhs(5)), strain, nel_strain)
   call mxCopyPtrToReal8(mxGetData(prhs(6)), dt, nel_dt)
   call mxCopyPtrToReal8(mxGetData(prhs(7)), tt, nel_tt)

   ! Determine a supported set of direct and shear components out of nel_strain
   call Estimate_Components(nel=nel_strain, num_dimensions=num_dimensions, num_shear=num_shear)
   if (.not. Check_Input_Dimensions(num_dimensions=num_dimensions, num_shear=num_shear)) then
      call Write_Error_And_Exit('mexFunction: ' // Formatval('nel_strain', nel_strain) &
         // ' can not be converted to a suitable set of direct and shear components')
   end if

   ! --- Assign to custom precision variables (Explicit type conversion)
   identifier = Uppercase(text=matname(1:setting_len_id))
   materialproperties = real(matprops, dp)
   totaltime = real(tt, dp)
   timeincrement = real(dt, dp)

   inp_stress = Import_Matrix(mat=real(stress, dp), num_dimensions=num_dimensions, num_shear=num_shear)
   tmp_dot_strain = Nonzero_Division(val=real(strain, dp), fac=timeincrement)
   inp_dot_strain = Import_Matrix(mat=tmp_dot_strain, num_dimensions=num_dimensions, num_shear=num_shear)
   inp_state = real(state, dp)

   ! --- Do calculation
   firstcall = .False.
   ! NOTE: Find a better method to determine the initial call of this routine for a simulation procedure
   if (totaltime <= timeincrement) then
      firstcall = .True.
   end if

   call Xmat_Initialize(xmat_obj=xmat_obj, &                         ! Assign variables internally for calculation
      solver_name=setting_solver_default, constitutive_model_name=identifier, &
      material_parameters=materialproperties(1:nel_props), calculateJacobian=.True., firstcall=firstcall)
   call xmat_obj%Import_State(nstates=nel_statev, &                  ! Prepare values for further calculation
      statevariables=inp_state(1:nel_statev), stress=inp_stress, dot_strain=inp_dot_strain, &
      totaltime=totaltime, timeincrement=timeincrement)
   call xmat_obj%Calculate_Results(exportstress=inp_stress, &        ! Do calculation and return values
      exportstates=inp_state, exportjacobian=jacobian_internal, exportdt=exportdt)

   ! --- Reassign to calling program variables (Implicit type conversion - might reducing precision of values)
   temp_mat = Export_Matrix(mat=inp_stress, num_dimensions=num_dimensions, num_shear=num_shear)
   stress = temp_mat
   state = inp_state
   jacobian = Export_Tensor(tens=jacobian_internal, num_dimensions=num_dimensions, num_shear=num_shear)

   ! Transpose jacobian before exporting to Matlab
   jacobian = transpose(jacobian)

   plhs(1) = mxCreateDoubleMatrix(nel_stress, 1, 0)
   plhs(2) = mxCreateDoubleMatrix(nel_statev, 1, 0)
   plhs(3) = mxCreateDoubleMatrix(nel_stress, nel_stress, 0)
   call mxCopyReal8ToPtr(stress, mxGetPr(plhs(1)), nel_stress)
   call mxCopyReal8ToPtr(state, mxGetPr(plhs(2)), nel_statev)
   call mxCopyReal8ToPtr(jacobian, mxGetPr(plhs(3)), nel_stress*nel_stress)
end subroutine mexFunction
