! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                              U M A T                               ! Entry point for Abaqus/Standard
! ------------------------------------------------------------------ !
subroutine UMAT( &
   ! ==================== user coding to define ==================== !
      stress, statev, ddsdde, sse, spd, scd, &
   ! ======================= and if necessary ====================== !
      rpl, ddsddt, drplde, drpldt, &
   ! ============= variables passed in for information ============= !
      stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
      ndi, nshr, ntens, nstatv, props, nprops, &
      coords, drot, &
   ! ================ variables that can be updated ================ !
      pnewdt, &
   ! ========== (more) variables passed in for information ========= !
      celent, dfgrd0, dfgrd1, &
      noel, npt, layer, kspt, jstep, kinc)
! ------------------------------------------------------------------ !
   use General_Settings, only: dp, setting_len_id, setting_solver_default, &
                               Uppercase, Check_Input_Dimensions, Write_Error_And_Exit
   use Debug, only: Formatval
   use Math_Operations, only: Nonzero_Division, Import_Matrix, Export_Matrix, Export_Tensor
   use Xmat_Class, only: Xmat, Xmat_Initialize
   use Custom_Utilities, only: Is_First_Call, Preprocess_State_Matrices, Postprocess_State_Matrices
   !
#ifdef ABQ_STD_CALLING
   include 'aba_param.inc'
#else
   implicit none

   real(dp) :: stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, &
               time, dtime, temp, dtemp, predef, dpred, props, coords, drot, pnewdt, celent, dfgrd0, dfgrd1
#endif
   !
   ! ============= variables passed in for information ============= !
   integer :: ndi                                                    ! Number of direct stress components at this point
   integer :: nshr                                                   ! Number of engineering sher stress components at this point
   integer :: ntens                                                  ! Size of the stress or strain component array (ndi+nshr)
   integer :: nstatv                                                 ! Number of solution-dependent state variables associated w. material
   integer :: nprops                                                 ! User-defined number of material constants associated w.this material
   dimension :: props(nprops)                                        ! User-specified array of material const. associated w. this material
   ! ==================== user coding to define ==================== !
   ! All those values passed as values at the beginning of the increment.
   ! They have to be updated in this routine and overwritten to be returned as end of increment values
   dimension :: stress(ntens)                                        ! stress tensor (as an array)
   dimension :: statev(nstatv)                                       ! Array containing solution-dependent state variables
   dimension :: ddsdde(ntens, ntens)                                 ! Jacobian matrix of the constitutive model `\frac{\partial\Delta\sigma}{\partial\Delta\epsilon}`
   !                                                                 ! where `\Delta\sigma` are the stress increments and `\Delta\epsilon` the strain increments
   !                                                                 ! Will be assumed symmetric unless otherwise indicated
   ! sse                                                             ! Specific elastic strain energy
   ! spd                                                             ! Specific plastic dissipation
   ! scd                                                             ! Specific "creep" dissipation
   ! ======= can be defined if necessary (thermal/pore fluid) ====== !
   ! rpl                                                             ! Volumetric heat generation per unit time caused
   !                                                                 ! by mechanical workings of material
   !                                                                 ! OR in geostatic stress procedure or coupled pore fluid analysis:
   !                                                                 ! rpl indicats if a cohesive element allows tangential pore fluid flow
   ! ============ can be defined if necessary (thermal) ============ !
   dimension :: ddsddt(ntens)                                        ! Variation of the stress increments with respect to temperature
   dimension :: drplde(ntens)                                        ! Variation of rpl with respect to strain increments
   ! drpldt                                                          ! Variation of rpl with respect to the temperature
   ! ============= variables passed in for information ============= !
   dimension :: stran(ntens)                                         ! Total strains at the beginning of the increment
   dimension :: dstran(ntens)                                        ! Array of strain increments
   dimension :: time(2)                                              ! time(1): value of step time at beginning of current increment/freq
   !                                                                 ! time(2): value of total time at the beginning of current increment
   ! dtime                                                           ! Time increment
   ! temp                                                            ! Temperature at the start of the increment
   ! dtemp                                                           ! Increment of temperature
   dimension :: predef(1)                                            ! Array of interpolated values of predefined field variables
   dimension :: dpred(1)                                             ! Array of increments of predefined field variables
   character(len=80) :: cmname                                       ! User-specified material name, always passed UPPERCASE.
   !                                                                 ! Don't use "ABQ_" as leading string
   dimension :: coords(3)                                            ! An array containing the coordinates of this point
   dimension :: drot(3, 3)                                           ! Rotation increment matrix
   ! ================ variables that can be updated ================ !
   ! pnewdt                                                          ! Ratio of suggested new time increment to time increment being used
   ! ========== (more) variables passed in for information ========= !
   ! celent                                                          ! Characteristic element length
   dimension :: dfgrd0(3, 3)                                         ! Array containing deformation gradient at the beginning of the incr.
   dimension :: dfgrd1(3, 3)                                         ! Array containing deformation gradient at the end of the increment
   integer :: noel                                                   ! Element number
   integer :: npt                                                    ! Integration point number
   integer :: layer                                                  ! Layer number (for composite shells and layered soils
   integer :: kspt                                                   ! Section point number within the current layer
   integer :: jstep(4)                                               ! jstep(1): step number
   !                                                                 ! jstep(2): procedure type key
   !                                                                 ! jstep(3): 1 if NLGEOM=YES for the current step, 0 otherwise
   !                                                                 ! jstep(4): 1 if current step is a linear perturbation procedure
   integer :: kinc                                                   ! Increment number
   ! --------------------------------------------------------------- !
   ! --- Variables with custom precision for the actual calculations
   real(dp), dimension(nprops) :: materialproperties
   real(dp), dimension(nstatv) :: inp_state
   real(dp), dimension(__matrix__) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(__tensor__) :: jacobian
   character(len=13), parameter :: intername = 'umat_abq_std '
   real(dp), dimension(ntens) :: tmp_stress, tmp_dot_strain
   real(dp), dimension(ntens, ntens) :: tmp_tens
   character(len=setting_len_id) :: identifier
   logical :: firstcall

   if (.not. Check_Input_Dimensions(num_dimensions=ndi, num_shear=nshr)) then
      call Write_Error_And_Exit('UMAT: This combination of ' // Formatval('ndi', ndi) // ' and ' // &
         Formatval('nshr', nshr) // ' is not supported')
   end if

#ifdef TOCHNOG_CALLING
   ! Tochnog material model is hardcoded at the beginning
   cmname = TOCHNOG_MATERIAL
#endif

   ! --- Assign to custom precision variables (consider element positioning)
   identifier = Uppercase(text=cmname(1:setting_len_id))

   materialproperties = real(props, dp)
   totaltime = real(time(2), dp)
   timeincrement = real(dtime, dp)

   inp_state = Preprocess_State_Matrices(identifier=identifier, intername=intername, nstates=nstatv, &
      states=real(statev, dp), rotation=real(drot, dp))

   tmp_stress = real(stress, dp)
   tmp_dot_strain = Nonzero_Division(val=real(dstran, dp), fac=timeincrement)

   ! For this interface (as required by Abaqus/Standard) engineering strain and shear components
   ! in the order (12, 13, 23) are expected for ndi = 3 and nshr = 3 (instead of (12, 23, 13))
   ! Adjust it for internal processing (and change it back before returning the results)
   if ((ndi == 3) .and. (nshr == 3)) then
      tmp_stress = [tmp_stress(1:4), tmp_stress(6), tmp_stress(5)]
      tmp_dot_strain = [tmp_dot_strain(1:4), tmp_dot_strain(6), tmp_dot_strain(5)]
   end if

   tmp_dot_strain(ndi+1:ntens) = 0.5_dp*tmp_dot_strain(ndi+1:ntens)

   inp_stress = Import_Matrix(mat=tmp_stress, num_dimensions=ndi, num_shear=nshr)
   inp_dot_strain = Import_Matrix(mat=tmp_dot_strain, num_dimensions=ndi, num_shear=nshr)

   firstcall = Is_First_Call(step=jstep(1), iteration=kinc)

   ! --- Do calculation
   call Xmat_Initialize(xmat_obj=xmat_obj, &                         ! Assign variables internally for calculation
      solver_name=setting_solver_default, constitutive_model_name=identifier, &
      material_parameters=materialproperties, calculateJacobian=.True., firstcall=firstcall)
   call xmat_obj%Import_State(nstates=nstatv, &                      ! Prepare values for further calculation
      statevariables=inp_state, stress=inp_stress, dot_strain=inp_dot_strain, &
      totaltime=totaltime, timeincrement=timeincrement)
   call xmat_obj%Calculate_Results(exportstress=inp_stress, &        ! Do calculation and return values
      exportstates=inp_state, exportjacobian=jacobian, exportdt=exportdt)

   ! --- Reassign to calling program variables (Implicit type conversion)
   tmp_stress = Export_Matrix(mat=inp_stress, num_dimensions=ndi, num_shear=nshr)

   if ((ndi == 3) .and. (nshr == 3)) then
      tmp_stress = [tmp_stress(1:4), tmp_stress(6), tmp_stress(5)]
   end if

   stress = tmp_stress

   inp_state = Postprocess_State_Matrices(identifier=identifier, intername=intername, nstates=nstatv, states=inp_state)
   statev = inp_state

   ! NOTE: Time increment updating not tested
   ! pnewdt = Nonzero_Division(val=exportdt, fac=timeincrement)

   tmp_tens = Export_Tensor(tens=jacobian, num_dimensions=ndi, num_shear=nshr)
   if ((ndi == 3) .and. (nshr == 3)) then
      tmp_tens(:, 5:6) = reshape([tmp_tens(:, 6), tmp_tens(:, 5)], [6, 2])
      tmp_tens(5:6, :) = transpose(reshape([tmp_tens(6, :), tmp_tens(5, :)], [6, 2]))
   end if

   ddsdde = tmp_tens
end subroutine UMAT
