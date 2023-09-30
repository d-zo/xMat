! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                        D o t _ V a l u e s                         !
! ------------------------------------------------------------------ !
subroutine Dot_Values( &
   ! ============= variables passed in for information ============= !
      materialname, nparams, materialparameters, nstatevar, statevariables, &
      ncomponents, oldstress, oldstrain, timeincrement, totaltime, &
   ! ==================== user coding to define ==================== !
      dotstress, dotstate, jacobian)
! ------------------------------------------------------------------ !
   use General_Settings, only: dp, setting_len_id, setting_solver_default, &
                               Estimate_Components, Uppercase, Check_Input_Dimensions, Write_Error_And_Exit
   use Debug, only: Formatval
   use Math_Operations, only: Nonzero_Division, Import_Matrix, Export_Matrix, Export_Tensor
   use Xmat_Class, only: Xmat, Xmat_Initialize
   implicit none
   !
   character(len=80), intent(in) :: materialname
   integer, intent(in) :: nparams, nstatevar
   real(dp), dimension(nparams), intent(in) :: materialparameters
   real(dp), dimension(nstatevar), intent(in) :: statevariables
   integer, intent(in) :: ncomponents
   real(dp), dimension(ncomponents), intent(in) :: oldstress, oldstrain
   real(dp), intent(in) :: timeincrement, totaltime
   real(dp), dimension(ncomponents), intent(out) :: dotstress
   real(dp), dimension(nstatevar), intent(out) :: dotstate
   real(dp), dimension(ncomponents, ncomponents), intent(out) :: jacobian
   ! --------------------------------------------------------------- !
   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(__tensor__) :: jacobian_internal
   real(dp), dimension(__matrix__) :: inp_stress, inp_dot_strain, dotstress_internal
   real(dp), dimension(ncomponents) :: tmp_dot_strain
   character(len=setting_len_id) :: identifier
   character(len=13), parameter :: intername = 'dot_values   '
   integer :: num_dimensions, num_shear
   logical :: firstcall

   ! Determine a supported set of direct and shear components out of ncomponents
   call Estimate_Components(nel=ncomponents, num_dimensions=num_dimensions, num_shear=num_shear)
   if (.not. Check_Input_Dimensions(num_dimensions=num_dimensions, num_shear=num_shear)) then
      call Write_Error_And_Exit('Dot_Values: ' // Formatval('ncomponents', ncomponents) &
         // ' can not be converted to a suitable set of direct and shear components')
   end if

   ! --- Assign to custom precision variables (consider element positioning)
   identifier = Uppercase(text=materialname(1:setting_len_id))

   inp_stress = Import_Matrix(mat=oldstress, num_dimensions=num_dimensions, num_shear=num_shear)
   tmp_dot_strain = Nonzero_Division(val=oldstrain, fac=timeincrement)
   inp_dot_strain = Import_Matrix(mat=tmp_dot_strain, num_dimensions=num_dimensions, num_shear=num_shear)

   firstcall = .False.
   ! NOTE: Find a better method to determine the initial call of this routine for a simulation procedure
   if (totaltime <= timeincrement) then
      firstcall = .True.
   end if

   call Xmat_Initialize(xmat_obj=xmat_obj, &                         ! Assign variables internally for calculation
      solver_name=setting_solver_default, constitutive_model_name=identifier, &
      material_parameters=materialparameters, calculate_jacobian=.True., firstcall=firstcall)
   call xmat_obj%Import_State(nstates=nstatevar, &                   ! Prepare values for further calculation
      statevariables=statevariables, stress=inp_stress, dot_strain=inp_dot_strain, &
      totaltime=totaltime, timeincrement=timeincrement)
   call xmat_obj%Dot_Results(dotstress=dotstress_internal, &         ! Do calculation and return values
      dotstates=dotstate, jacobian=jacobian_internal)

   ! --- Reassign to calling program variables
   dotstress = Export_Matrix(mat=dotstress_internal, num_dimensions=num_dimensions, num_shear=num_shear)
   jacobian = Export_Tensor(tens=jacobian_internal, num_dimensions=num_dimensions, num_shear=num_shear)
end subroutine Dot_Values
