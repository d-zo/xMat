! ==================================================================================================================== !
module General_Settings
   implicit none

   integer, parameter :: dp = selected_real_kind(15)                 ! Desired precision for calculation used everywhere in this routine.
   !                                                                 ! Double precision (15) by default, higher values might be possible.
   real(dp), parameter :: setting_epsilon = 500.0_dp*epsilon(1.0_dp) ! Very small number (only 3 orders away from machine precision). Used
   !                                                                 ! to check if numbers are so small to be regarded as zero

   ! --- Constitutive routines ---
   integer, parameter :: setting_len_id = 9                          ! Length and name of identifiers for all available constitutive models
   ! NOTE: Although the interfaces allow for case insensitive identifiers, the setting_id_... variables have to be all uppercase
   character(len=setting_len_id), parameter :: setting_id_elasticity     = 'ELAS-0000'
   character(len=setting_len_id), parameter :: setting_id_hypo_wu92      = 'HYPO-WU92'
   character(len=setting_len_id), parameter :: setting_id_hypo_vw96      = 'HYPO-VW96'
   character(len=setting_len_id), parameter :: setting_id_viscohypo_ni03 = 'VIHY-NI03'
   character(len=setting_len_id), parameter :: setting_id_barodesy_ko15  = 'BARO-KO15'
   character(len=setting_len_id), parameter :: setting_id_barodesy_sc18  = 'BARO-SC18'
   character(len=setting_len_id), parameter :: setting_id_barodesy_ko21  = 'BARO-KO21'
   character(len=setting_len_id), parameter :: setting_id_test_dgl       = 'TEST-DGLX'
   !
   ! Tochnog does not support selecting a material model by passing its name. The used material modes must be hardcoded here
#ifdef TOCHNOG_CALLING
#define TOCHNOG_MATERIAL 'HYPO-VW96'
#endif

   ! --- Jacobian ---
   logical, parameter :: setting_numerical_jacobian = .False.        ! Approximate the jacobian numerically, even if the constitutive model
   !                                                                 ! could provide a jacobian. If set to .False., the jacobian of the
   !                                                                 ! constitutive model will be used. In case the constitutive model does
   !                                                                 ! not provide a jacobian, the numerical approximation will be used as
   !                                                                 ! as if this value were set to .True.
   real(dp), parameter :: setting_numerical_jacobian_disturbance = 0.00001_dp
   !                                                                 ! The disturbance used for the numerical approximation of the jacobian
   !                                                                 ! Should not be too small (the number of significant digits of the
   !                                                                 ! jacobian is less than the number of orders between this and epsilon)

   ! --- Hypoplasticity ---
   real(dp), parameter :: setting_very_small_stress = -0.01_dp       ! Upper limit for stresses (to issue a "close to tension"-warning)
   real(dp), parameter :: setting_min_youngs_modulus = 5.0_dp        ! Lower limit for (approximated) Young's modulus. A high stiffness
   !                                                                 ! results in a bigger difference when transitioning between elastic
   !                                                                 ! and hypoplastic response but yields larger steps in tensile sections
   logical, parameter :: setting_hypo_consistent_f_d = .True.        ! Use a correction term `\bar{f}_d` for a consistent lower bound of `f_d`
   logical, parameter :: setting_hypo_increased_stiffness = .False.  ! Modifications to calculation of `\mathcal{L}`, `\mathbf{N}` and therefore `\overset{\circ}{\mathbf{T}}`

   ! --- Solver ---
   character(len=7), parameter :: setting_solver_default = 'Eul-Exp' ! Specify the default solver name. Currently available names are
   !                                                                 ! 'Eul-Exp', 'Richard', 'RK-S 23', 'RK-BS23', 'RK-F 45', 'RK-DP45' and
   !                                                                 ! 'RK-CK45'
   integer, parameter :: setting_max_integration_steps = 50000       ! This is the maximum number of steps allowed for a given time
   !                                                                 ! increment.
   logical, parameter :: setting_restrict_initial_substep = .True.   ! Uses setting_initial_substep_scale for the first step (or all steps
   !                                                                 ! if integration method is non-adaptive) if True. Otherwise the given
   !                                                                 ! dt of the calling program is used (resulting in a single step for
   !                                                                 ! non-adaptive methods)
   real(dp), parameter :: setting_initial_substep_scale = 0.00001_dp ! If enforced with setting_enforce_initial_substep_dt the minimum of
   !                                                                 ! this and the given dt will be used for the (first) time increment
   real(dp), parameter :: setting_stepsize_scaling_safety = 0.9      ! Safety factor applied to the estimation of new time increments
   !                                                                 ! during the integration loop (see also Press, 1997, p. 712)
   real(dp), parameter :: setting_max_rel_error = 0.001_dp           ! Maximum relative error for adaptive integration methods. The error
   !                                                                 ! is scaled componentwise first and then the maximum absolute error
   !                                                                 ! has to be smaller than this value (see also Integrate subroutine in
   !                                                                 ! the solver module)
   integer, parameter :: setting_max_integration_refinements = 20    ! When the error is too large, an integration step is repeated with a
   !                                                                 ! smaller time increment (and checked again). Stop if those checks
   !                                                                 ! fail repeatedly without any accepted

   ! --- Internal memory reservation for Matlab (mexFunction) and Plaxis (USER_MOD). Might need to be adjusted for new constitutive models
   integer, parameter :: setting_num_statevariables = 20             ! Maximum number of state variables for each constitutive model
   integer, parameter :: setting_max_mat_params = 16                 ! Maximum number of provided material parameters
   integer, parameter :: setting_max_internal_states = &             ! Internal states include stresses (+jac) and statevariables (+jac)
                                                       (__nelmat__+1)*__nelmat__ + (__nelmat__+1)*setting_num_statevariables

   ! --- Logs ---
   character(len=6), parameter :: identifier_log_start = 'xMat: '    ! Identifier written at the beginning of each output by write_out
   !                                                                 ! (is only used by functions within this module)
   !
   ! --- Global variables ---
   ! ATTENTION: These variables should only be changed by functions of this module (Don't change them unless you know what you're doing)
   integer :: global_num_direct_components = 3                       ! Default to three until actual values are read (2 or 3)
   integer :: global_num_shear_components = 3                        ! Default to three until actual values are read (0 to 3, depending
                                                                     ! on element type)

   private
   public dp, setting_epsilon, setting_len_id, setting_id_elasticity, setting_id_hypo_wu92, setting_id_hypo_vw96, &
          setting_id_viscohypo_ni03, setting_id_barodesy_ko15, setting_id_barodesy_sc18, setting_id_barodesy_ko21, &
          setting_id_test_dgl, &
          setting_hypo_consistent_f_d, setting_hypo_increased_stiffness, setting_very_small_stress, &
          setting_min_youngs_modulus, setting_solver_default, setting_max_integration_steps, &
          setting_restrict_initial_substep, setting_initial_substep_scale, &
          setting_stepsize_scaling_safety, setting_max_rel_error, setting_max_integration_refinements, &
          setting_numerical_jacobian, setting_numerical_jacobian_disturbance, &
          setting_num_statevariables, setting_max_mat_params, setting_max_internal_states, &
          global_num_direct_components, global_num_shear_components, &
          Estimate_Components, Check_Input_Dimensions, Is_Nan, C_To_Fortran_String, Uppercase, &
          Write_Message, Write_Warning, Write_Error_And_Exit


   contains


#addfile subroutine Estimate_Components


#addfile function Check_Input_Dimensions


#addfile function Is_Nan


#addfile function C_To_Fortran_String


#addfile function Uppercase


#addfile subroutine Write_Message


#addfile subroutine Write_Warning


#addfile subroutine Write_Error_And_Exit
end module General_Settings
