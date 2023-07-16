! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                          U S E R _ M O D                           ! Entry point for Plaxis
! ------------------------------------------------------------------ !
subroutine USER_MOD( &
   ! ============= variables passed in for information ============= !
      IDTask, iMod, IsUndr, iStep, iTer, Iel, Jnt, X, Y, Z, &
      Time0, dTime, Props, Sig0, Swp0, StVar0, dEps, &
   ! ==================== user coding to define ==================== !
      D, Bulk_W, Sig, Swp, StVar, &
      ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang, &
   ! ========== variables passed in for debug information ========== !
      iPrjDir, iPrjLen, &
   ! ==================== user coding to define ==================== !
      iAbort)
! ------------------------------------------------------------------ !
   use General_Settings, only: dp, setting_len_id, setting_max_mat_params, setting_solver_default, &
                               Check_Input_Dimensions, Write_Error_And_Exit
   use Math_Operations, only: Nonzero_Division, Import_Matrix, Export_Matrix, Export_Tensor
   use PlaxisInformationPool, only: ModelList, numParameters, numStateVariables
   use Xmat_Class, only: Xmat, Xmat_Initialize
   use Custom_Utilities, only: Is_First_Call, Preprocess_State_Matrices, Postprocess_State_Matrices
   implicit none
   !
   ! ============= variables passed in for information ============= !
   integer :: IDTask                                                 ! Identification of task where
   !                                                                 !   1: Initialise state variables
   !                                                                 !   2: Calculate constitutive stresses
   !                                                                 !   3: Create effective material stiffness matrix
   !                                                                 !   4: Return number of state variables
   !                                                                 !   5: Return matrix attributes
   !                                                                 !   6: Create elastic material stiffness matrix
   integer :: iMod                                                   ! User defined soil model number as defined in Modellist from
   !                                                                 ! PlaxisInformationPool
   integer :: nStat                                                  ! Number of state variables
   integer :: IsUndr                                                 ! 0: Drained condition; 1: Undrained condition
   integer :: iStep                                                  ! Current calculation step number
   integer :: iTer                                                   ! Current iteration number
   integer :: Iel                                                    ! Current element number
   integer :: Jnt                                                    ! Current local stress point number (1 to 3 for 6-noded elements)
   double precision :: X, Y, Z                                       ! Global coordinates of current stress point
   double precision :: Time0                                         ! Time at start of current step
   double precision :: dTime                                         ! Time increment of current step
   !
   ! --- Variables regarding current stress point ---
   double precision, dimension(50) :: Props                          ! Array with user-defined model parameters
   double precision, dimension(20) :: Sig0                           ! Array with
   !                                                                 !  (1-6) effective stress components at start of current step
   !                                                                 !  (7-13) `p_{steady}`, `\sum Mstage^0`, `\sum Mstage`, `Sat`, `Sat^0`, `Suc`, `Suc^0`
   !                                                                 !  (14-20) `\sum Msf^0`, `\sum Msf`, `X_j`, SatRes, Temp, UnfrozenW, 0
   double precision :: Swp0                                          ! Previous excess pore pressure
   ! IMPORTANT: For IDTask = 1, StVar0 is an input and output parameter
   double precision, dimension(nStat) :: StVar0                      ! Array with values of state variables
   double precision, dimension(12) :: dEps                           ! Array with strain increments (1-6) and
   !                                                                 ! initial strains (7-12) at start of current step
   ! ==================== user coding to define ==================== !
   double precision, dimension(6, 6) :: D                            ! Effective material stiffness matrix
   double precision :: Bulk_W                                        ! Bulk modulus of water (for undrained calculations and consolidation)
   double precision, dimension(6) :: Sig                             ! Array with resulting constitutive stresses
   double precision :: Swp                                           ! Resulting excess pore pressure
   double precision, dimension(nStat) :: StVar                       ! Array with resulting values of state variables
   !
   ! --- General indicators
   integer :: ipl                                                    ! Plasticity indicator:
   !                                                                 !   0: no plasticity
   !                                                                 !   1: Mohr-Coulomb (failure) point
   !                                                                 !   2: Tension cut-off point
   !                                                                 !   3: Cap hardening point
   !                                                                 !   4: Cap friction point
   !                                                                 !   5: Friction hardening point
   integer :: NonSym                                                 ! Stiffness matrix 0: is symmetric; 1: is not symmetric
   integer :: iStrsDep                                               ! Stiffness matrix 0: is not stress-dependent; 1: is stress dependent
   integer :: iTimeDep                                               ! Stiffness matrix 0: is not time-dependent; 1: is time-dependent
   integer :: iTang                                                  ! Stiffness matrix
   !                                                                 !   0: is not a tangent stiffness
   !                                                                 !   1: is a tangent stiffness matrix
   !                                                                 ! to be used in a full Newton-Raphson iteration process
   ! ========== variables passed in for debug information ========== !
   integer :: iPrjDir                                                ! Project directory
   integer :: iPrjLen                                                ! Length of project directory name
   ! ==================== user coding to define ==================== !
   integer :: iAbort                                                 ! Forcing calculation to stop if set to 1
   ! --------------------------------------------------------------- !
   !implicit double precision (a-h, o-z)

   ! --- Variables with custom precision for the actual calculations
   real(dp), dimension(setting_max_mat_params) :: materialproperties
   real(dp), dimension(nStat) :: inp_state
   real(dp), dimension(__matrix__) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(__tensor__) :: jacobian
   character(len=setting_len_id) :: identifier
   character(len=13), parameter :: intername = 'user_mod     '
   real(dp), dimension(6) :: tmp_dot_strain                          ! Assume Plaxis input to be vec6
   real(dp), parameter, dimension(3, 3) :: rotation = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
   integer :: num_dimensions, num_shear, nparams
   logical :: firstcall

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: User_Mod             ! Export function name in dll
#endif

   if ((iMod >= 1) .and. (iMod <= len(ModelList))) then
      identifier = ModelList(iMod)
   else
      iAbort = 1
      call Write_Error_And_Exit('USER_MOD: Requested soil model not present')
   end if

   if (IDTask == 4) then
      nStat = numStateVariables(iMod)
      return
   end if

   if (IDTask == 5) then
      NonSym   = 1
      iStrsDep = 1
      iTang    = 1
      iTimeDep = 0
      return
   end if

   ! Default output variables
   iAbort = 0
   ipl = 0

   ! Assuming that Plaxis is always submitting vectors with 3 dimensions and 3 shear components (otherwise change definitions above)
   num_dimensions = 3
   num_shear = 3
   if (.not. Check_Input_Dimensions(num_dimensions=num_dimensions, num_shear=num_shear)) then
      call Write_Error_And_Exit('USER_MOD: This combination of direct/shear components is not supported')
   end if

   ! Get number of material parameters from chosen constitutive model and data saved in PlaxisInformationPool
   nparams = numParameters(iMod)

   ! --- Assign to custom precision variables (Explicit type conversion)
   ! For this interface (as required by Plaxis) engineering strain is required. Adjust it for internal processing
   totaltime = real(Time0, dp)
   timeincrement = real(dTime, dp)

   inp_stress = Import_Matrix(mat=real(Sig0(1:6), dp), num_dimensions=num_dimensions, num_shear=num_shear)

   tmp_dot_strain = Nonzero_Division(val=[real(dEps(1:3), dp), 0.5_dp*real(dEps(4:6), dp)], fac=timeincrement)
   inp_dot_strain = Import_Matrix(mat=tmp_dot_strain, num_dimensions=num_dimensions, num_shear=num_shear)

   materialproperties(1:nparams) = real(Props(1:nparams), dp)

   if (IDTask == 1) then
      inp_state = real(Props((nparams+1):(nparams+nStat)), dp)
   else
      inp_state = real(StVar0, dp)
   end if
   inp_state = Preprocess_State_Matrices(identifier=identifier, intername=intername, nstates=nStat, &
      states=inp_state, rotation=rotation)

   firstcall = Is_First_Call(step=iStep, iteration=iTer)

   ! --- Do calculation
   call Xmat_Initialize(xmat_obj=xmat_obj, &                         ! Assign variables internally for calculation
      solver_name=setting_solver_default, constitutive_model_name=identifier, &
      material_parameters=materialproperties(1:nparams), calculateJacobian=.True., firstcall=firstcall)
   call xmat_obj%Import_State(nstates=nStat, &                       ! Prepare values for further calculation
      statevariables=inp_state, stress=inp_stress, dot_strain=inp_dot_strain, &
      totaltime=totaltime, timeincrement=timeincrement)
   call xmat_obj%Calculate_Results(exportstress=inp_stress, &        ! Do calculation and return values
      exportstates=inp_state, exportjacobian=jacobian, exportdt=exportdt)

   inp_state = Postprocess_State_Matrices(identifier=identifier, intername=intername, nstates=nStat, states=inp_state)

   ! --- Reassign custom precision variables to corresponding output variables (Implicit type conversion)
   if (IDTask == 1) then
      StVar0 = inp_state
   end if

   if (IDTask == 2) then
      Sig = Export_Matrix(mat=inp_stress, num_dimensions=num_dimensions, num_shear=num_shear)
      StVar = inp_state
      Swp = Swp0
   end if

   if ((IDTask == 3) .or. (IDTask == 6)) then
      D = Export_Tensor(tens=jacobian, num_dimensions=num_dimensions, num_shear=num_shear)
   end if
end subroutine USER_MOD
