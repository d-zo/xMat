! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                             V U M A T                              ! Entry point for Abaqus/Explicit
! ------------------------------------------------------------------ !
subroutine VUMAT( &
   ! ========= (scalar) variables passed in for information ======== !
      nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &       ! integer
      stepTime, totalTime, dt, &                                     ! real
   ! ========= (other) variables passed in for information ========= !
      cmname, &                                                      ! character
      coordMp, charLength, props, density, strainInc, relSpinInc, &
      tempOld, stretchOld, defgradOld, fieldOld, stressOld, stateOld, enerInternOld, enerInelasOld, &
      tempNew, stretchNew, defgradNew, fieldNew, &
   ! ==================== user coding to define ==================== !
      stressNew, stateNew, &
   ! ================ variables that can be updated ================ !
      enerInternNew, enerInelasNew)
! ------------------------------------------------------------------ !
   use General_Settings, only: dp, setting_len_id, setting_solver_default, &
                               Uppercase, Check_Input_Dimensions, Write_Error_And_Exit
   use Debug, only: Formatval
   use Math_Operations, only: Nonzero_Division, Import_Matrix, Export_Matrix
   use Xmat_Class, only: Xmat, Xmat_Initialize
   !
#ifdef ABQ_EXP_CALLING
   ! ATTENTION: Using Abaqus/2018 the first seven comment lines in vaba_param_sp.inc and vaba_param_dp.inc
   !            have to be adjusted (! instead of c as comment start) or the include will produce syntax errors.
   !            The files are typically located in Abaqus/SimulationServices/V6R2018x/SMAUsubs/PublicInterfaces
   include 'vaba_param.inc'
#else
   implicit none

   real(dp) :: stepTime, totalTime, dt, coordMp, charLength, props, density, strainInc, relSpinInc, &
      tempOld, stretchOld, defgradOld, fieldOld, stressOld, stateOld, enerInternOld, enerInelasOld, &
      tempNew, stretchNew, defgradNew, fieldNew, stressNew, stateNew, enerInternNew, enerInelasNew
#endif
   !
   ! ========= (scalar) variables passed in for information ======== !
   integer, intent(in) :: nblock                                     ! Number of material points to be processed in this call
   integer, intent(in) :: ndir                                       ! Number of direct components in a symmetric tensor
   integer, intent(in) :: nshr                                       ! Number of indirect components in a symmetric tensor
   integer, intent(in) :: nstatev                                    ! Number of user-defined state variables associated with material type
   integer, intent(in) :: nfieldv                                    ! Number of user-defined external variable fields
   integer, intent(in) :: nprops                                     ! User-specified number of user-defined material properties
   integer, intent(in) :: lanneal                                    ! Flag indication whether the routine is being called
   !                                                                 ! during an annealing process (0=normal mechanics)
   ! stepTime                                                        ! Value of time since the step began
   ! totalTime                                                       ! Value of total time. The time at the beginning of the step
   !                                                                 ! is given by totalTime-stepTime
   ! dt                                                              ! Time increment size
   ! ========= (other) variables passed in for information ========= !
   character(len=80) :: cmname                                       ! User-specified material name, always passed UPPERCASE.
   !                                                                 ! Don't use "ABQ_" as leading string
   dimension :: coordMp(nblock, *)                                   ! Material point coordinates
   dimension :: charLength(nblock)                                   ! Characteristic element length. Default value is typical length
   !                                                                 ! of a line across an element
   dimension :: props(nprops)                                        ! User-supplied material properties
   dimension :: density(nblock)                                      ! Current density at the material points in the midstep configuration
   !                                                                 ! (may be inaccurate when volumentric strain increment is very small,
   !                                                                 ! not affected by mass scaling)
   dimension :: strainInc(nblock, ndir+nshr)                         ! Strain increment tensor at each material point
   dimension :: relSpinInc(nblock, nshr)                             ! Incremental relative rotation vector at each material point
   !                                                                 ! defined in the corotational system
   ! --- Values at the beginning of the increment (all with suffix Old) ---
   dimension :: tempOld(nblock)                                      ! Temperatures at each material point
   dimension :: stretchOld(nblock, ndir+nshr)                        ! Stress tensor U at each material point defined from the
   !                                                                 ! polar decomposition of the deformation gradient by `\mathbf{F}=\mathbf{R}\cdot \mathbf{U}`
   dimension :: defGradOld(nblock, ndir+nshr+nshr)                   ! Deformation gradient tensor at each material point (stored as 9 (3D)
   !                                                                 ! or 5 (2D) element vector)
   dimension :: fieldOld(nblock, nfieldv)                            ! Values of the user-defined field variables at each material point
   dimension :: stressOld(nblock, ndir+nshr)                         ! Stress tensor at each material point
   dimension :: stateOld(nblock, nstatev)                            ! State variables at each material point
   dimension :: enerInternOld(nblock)                                ! Internal energy per unit mass at each material point
   dimension :: enerInelasOld(nblock)                                ! Dissipated inelastic energy per unit mass at each material point
   ! Values at the end of the increment (all with suffix New)
   dimension :: tempNew(nblock)                                      ! Temperatures at each material point
   dimension :: stretchNew(nblock, ndir+nshr)                        ! Stress tensor `\mathbf{T}` at each material point defined from the
   !                                                                 ! polar decomposition of the deformation gradient by `\mathbf{F}=\mathbf{R}\cdot \mathbf{U}`
   dimension :: defgradNew(nblock, ndir+nshr+nshr)                   ! Deformation gradient tensor at each material point (stored as 9 (3D)
   !                                                                 ! or 5 (2D) element vector)
   dimension :: fieldNew(nblock, nfieldv)                            ! Values of the user-defined field variables at each material point
   ! ==================== user coding to define ==================== !
   dimension :: stressNew(nblock, ndir+nshr)                         ! Stress tensor at each material point
   dimension :: stateNew(nblock, nstatev)                            ! State variables at each material point
   ! ================ variables that can be updated ================ !
   dimension :: enerInternNew(nblock)                                ! Internal energy per unit mass at each material point
   dimension :: enerInelasNew(nblock)                                ! Dissipated inelastic energy per unit mass at each material point
   ! --------------------------------------------------------------- !
   !! --- Variables with custom precision for the actual calculations
   real(dp), dimension(nprops) :: materialproperties
   real(dp), dimension(nstatev) :: inp_state
   real(dp), dimension(__matrix__) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, steptime_internal, totaltime_internal, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(__tensor__) :: jacobian
   character(len=13), parameter :: intername = 'vumat_abq_exp'
   real(dp), dimension(ndir+nshr) :: temp_mat, tmp_dot_strain
   character(len=setting_len_id) :: identifier
   logical :: firstcall
   integer :: idx_block

   if (.not. Check_Input_Dimensions(num_dimensions=ndir, num_shear=nshr)) then
      call Write_Error_And_Exit('VUMAT: This combination of ' // Formatval('ndi', ndir) // ' and ' // &
         Formatval('nshr', nshr) // ' is not supported')
   end if

   ! --- Assign to custom precision variables
   identifier = Uppercase(text=cmname(1:setting_len_id))
   materialproperties = real(props, dp)
   steptime_internal = real(stepTime, dp)
   totaltime_internal = real(totalTime, dp)
   timeincrement = real(dt, dp)

   ! --- Do calculation
   firstcall = .False.
   ! NOTE: Find a better method to determine the initial call of this routine for a simulation procedure
   if (totaltime_internal <= steptime_internal) then
      firstcall = .True.
   end if

   call Xmat_Initialize(xmat_obj=xmat_obj, &                         ! Assign variables internally for calculation
      solver_name=setting_solver_default, constitutive_model_name=identifier, &
      material_parameters=materialproperties, calculate_jacobian=.False., firstcall=firstcall)

   do idx_block = 1, nblock                                          ! Assign to custom precision variables
      inp_stress = Import_Matrix(mat=real(stressOld(idx_block, :), dp), num_dimensions=ndir, num_shear=nshr)
      tmp_dot_strain = Nonzero_Division(val=real(strainInc(idx_block, :), dp), fac=timeincrement)
      inp_dot_strain = Import_Matrix(mat=tmp_dot_strain, num_dimensions=ndir, num_shear=nshr)

      inp_state  = real(stateOld(idx_block, :), dp)

      call xmat_obj%Import_State(nstates=nstatev, &                  ! Prepare values for further calculation
         statevariables=inp_state, stress=inp_stress, dot_strain=inp_dot_strain, &
         totaltime=totaltime_internal, timeincrement=timeincrement)
      call xmat_obj%Calculate_Results(exportstress=inp_stress, &     ! Do calculation and return values
         exportstates=inp_state, exportjacobian=jacobian, exportdt=exportdt)

      ! --- Reassign to calling program variables (output types should not be stated explicitly)
      temp_mat = Export_Matrix(mat=inp_stress, num_dimensions=ndir, num_shear=nshr)
      stressNew(idx_block, :) = temp_mat
      stateNew(idx_block, :) = inp_state
   end do
end subroutine VUMAT
