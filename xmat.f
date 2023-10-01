!                                                     xMat constitutive model container
!  External program                                         D. Zobel (2023-09-30)                                         External program
!  _||_                                                                                                                                /\
!  \  /                                            released under GPL version 3 or greater                                            /  \
!   \/                                              see also https://github.com/d-zo/xMat                                              ||
!  +-------------------------------------------------------------------------------------------------------------------------------------+
!  |   I n t e r f a c e                                                                                                                 |
!  |=====================================================================================================================================|
!  |                                                                                                                                     |
!  | UMAT/VUMAT for Abaqus (6.14-2021), mexFunction for Matlab (2011b-2020b), USER_MOD for Plaxis (2015-2019) and                        |
!  | xmat_console as a simple Fortran test interface (for usage of the individual interfaces see below)                                  |
!  +-------------------------------------------------------------------------------------------------------------------------------------+
!                    ^                            ||                  ||                  ||                                        ^
!                    |                      Xmat_Initialize      Import_State      Calculate_Results                                |
!                    |                           .||.                .||.                .||.                                   export...
!                    |                            \/                  \/                  \/                                        |
!  +-------------------------------------+     +-----------------------------------------------------------------------------------------+
!  |   C u s t o m _ U t i l i t i e s   |     |   X m a t                                                                               |
!  |=====================================|     |=========================================================================================|
!  |                                     |     |                                                                                         |
!  | Auxiliary interface functions       |     | Store/load variables and control program flow                                           |
!  +-------------------------------------+     +-----------------------------------------------------------------------------------------+
!                                                   |               |                         :                                  ^
!  +-------------------------------------+          |    Select_Constitutive_Model        Integrate                              |
!  |   G e n e r a l _ S e t t i n g s   |          |               V                         :                                  |
!  |=====================================|          |          +---------------------------------------------------------+       |
!  |                                     |          |          |   C o n s t i t u t i v e _ M o d e l                   |       |
!  | 'setting_'-parameters encouraged    |          |          |=========================================================|       |
!  |            to be modified at will   |     Select_Solver   |                                                         |   new_state
!  |                                     |          |          | e.g. Hypoplasticity, Barodesy, Elasticity, ...          |       |
!  |  'global_'-variables are assigned   |          |          | 'param_'- and 'derived_'-values used as class variables |       |
!  |            later on and should not  |          |          +---------------------------------------------------------+       |
!  |            be edited here directly  |          |                                         :       ^            |             |
!  +-------------------------------------+          |                                         :    cur_state    dot_state        |
!                                                   V                                         V       |            V             |
!  +-------------------------------------+     +-----------------------------------------------------------------------------------------+
!  |   M a t h _ O p e r a t i o n s     |     |   S o l v e r                                                                           |
!  |=====================================|     |=========================================================================================|
!  |                                     |     |                                                                                         |
!  | 'const_'-parameters for mathe-      |     | Euler-Explicit, Runge-Kutta 2/3, ...                                                    |
!  |          matical expressions        |     |                                                                                         |
!  +-------------------------------------+     +-----------------------------------------------------------------------------------------+
!
! -----------------------
! -   O v e r v i e w   -
! -----------------------
! This file provides a container for different constitutive models and solvers. Similar functionalities are bundled in classes/modules.
! The choice, which constitutive model/solver configuration is used is controlled by (data from) the calling program and modifiable
! parameters in "General_Settings" module.
!
! This routine works with symmetric input matrices only. The default order of the components is as follows but certain interfaces use a
! different order of components (e.g. for UMAT. See respective interface definitions in this file):
!  - symmetric 3x3 matrices are stored as 6x1 vectors in Voigt notation.
!      The order of index pairs is `(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 3)`, i.e. a matrix `\mathbf{M}` becomes a vec6 `v^T=\left(M_{11}, M_{22}, M_{33}, M_{12}, M_{23}, M_{13}\right)`.
!  - symmetric 3x3x3x3 tensors are stored as 6x6 matrices, with tensorial index pairs `(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 3)` in both dimensions,
!      such that the (3, 4) entry of the 6x6 matrix represents the (3, 3, 1, 2) entry of the full tensor (and its symmetric entries).
!
! ---------------------------------
! -   H o w   t o   m o d i f y   -
! ---------------------------------
! All parameters encouraged to modify are stored in "General_Settings" module and called "setting_...". They are stored at the beginning of
! this file right after this introductory text.
!
! If actual program code should be extended (or modified due to erroneous behaviour), it is recommended to work with the official source
! files and not within this file ! They can be found at https://github.com/d-zo/xMat
! This repository contains all functions in organized and separate files and a tool to bundle them to a new single (big) file like this.
!
! ---------------------------------
! -   H o w   t o   e x t e n d   -
! ---------------------------------
! Since this library separates solver, constitutive routines and interfaces into modules (classes), an extension is rather straightforward:
!
! Adding a new constitutive model
!    1) Add a new module/class which extends Constitutive_Model_Baseclass and defines "Initialize()" and a "Get_Dot_State()". The interface
!       and typical usage of both can be seen from the existing implementations. If any other functions are used specifically to this
!       constitutive model, they can be added to this module. To use functions among multiple constitutive models, consider putting them in
!       the Constitutive_Model_Baseclass itself.
!    2) Select a suitable name and create a new "setting_id_..." variable in the "General_Settings" module.
!    3) Adjust "Select_Constitutive_Model()" subroutine, i.e. add "use" statement for class (and inclusion of files) as well as an
!        if-branch with an allocation command (if the material name is matched).
!    4) Add information about the new constitutive model to the values in the PlaxisInformationPool module.
!
! Adding a new solver
!    1) Add a new module/class which extends Solver_Baseclass and defines "Integration_Algorithm()". The interface and typical usage of
!       both can be seen from the existing implementations.
!    2) Adjust "Select_Solver()" subroutine, i.e. add "use" statement for class (and inclusion of files) as well as an if-branch with an
!       allocation command (if the solver name is matched). Currently one solver has to be chosen manually, so modify the variable
!       "setting_solver_default" in "General_Settings" module if needed.
!
! Adding a new interface
!    1) Add subroutine at the end of this file which uses "Xmat" and "Xmat_Initialize()" from the "Xmat_Class" module.
!    2) Instantiate an allocatable variable of "class(Xmat)".
!    3) Use the "Import_State()" and "Calculate_Results()" member function calls with this variable to do the actual assignments and
!       calculations. A typical usage can be seen from the existing implementations.
!    4) Make sure to convert input and output variables accordingly (i.e. dimensions and precision of expected variables).
!
! Adding custom interface utility functions
!    1) Add new functions/subroutines to "Custom_Utilities" module and import it in the required interface(s).
!
! -----------------------------------------------------------
! -   H o w   t o   c o m p i l e   ( a s   o b j e c t )   -
! -----------------------------------------------------------
! Standalone compilation (in Linux environment -- for Windows replace argument initialiser "-" with "/"):
!       gfortran -ffree-form -std=f2003 -cpp -o xmat.o -c xmat.f        ! for optimization add e.g.: -O3 -march=native -flto -Ofast
!       ifort -free -fpp -stand f03 -o xmat.o -c xmat.f                 ! for optimization add e.g.: -O3 -xHost -ipo -fast
!
! Shared library with C interface (in Linux environment):
!       gfortran -ffree-form -std=f2003 -cpp -shared -fPIC -DCINTER -o xmat.so xmat.f
!
! MATLAB (first command tested for Matlab 2015a and 2018a in Linux environment with GNU Fortran Compiler,
!    second command tested for Matlab 2019b in Windows environment with Intel Fortran Compiler):
!       mex FFLAGS='-ffree-form $FFLAGS' -largeArrayDims -DMATLAB_CALLING xmat.f
!       mex 'COMPFLAGS="$COMPFLAGS /free"' -largeArrayDims -DMATLAB_CALLING xmat.f
!
! Abaqus 6.14-2021 (in Linux environment, for Windows replace argument initialiser "-" with "/"):
!    First, the abaqus_v6.env file has to be modified: Add '-free' and '-DABAQUS_CALLING' as arguments to compile_fortran). Either edit
!    the system file (global change) or a local copy of abaqus_v6.env (local effect only). Abaqus also supports precompiling the routine:
!       abaqus make library=xmat.f
!
! Plaxis 2016 (Plaxis expects a dynamic linked library (.dll). In Windows open Intel Compiler command prompt or command window with
!    necessary library/PATH directives for compilation (architecture (32/64bit) must fit to Plaxis installation). Compile with the
!    following command and move/copy xmat64.dll to udsm-directory within Plaxis installation folder (2D and 3D) afterwards.
!       ifort /free /fpp /DPLAXIS_DLL /dll xmat.f /recursive /MT /exe:xmat64.dll
!
! If Plaxis functions or the C interface should be used directly and not as a library, pass the flag '-DNOBIB' at compilation
! to prevent conflicts with the '!DEC$' compiler directives.
!
! -------------------------------------
! -   T r o u b l e s h o o t i n g   -
! -------------------------------------
! If user routine compilation fails with "Illegal character in statement label field":
!    -> compile with free format, i.e. add '-free' to compile command (e.g. in Abaqus environment file).
! If user routine compilation fails with "Bad # preprocessor line":
!    -> compile with active preprocessing, i.e. add '-fpp' (gfortran) or '-cpp' (ifort) to compile command.
! If linking of user routine fails with "Undefined reference to ...":
!    -> If the reference is a xMat function, possibly the respective preprocessor macro is not defined ('-D<Macro>')
!       or the function has a '!DEC$' directive but is directly called (missing '-DNOBIB').
!    -> If the reference is an external function, possibly path ('-L<path>') or library switch ('-l<name>') are missing
! If an Abaqus simulation fails with "Abqaqus Error: Problem during linking".
!    -> Possibly macros for more interfaces than needed are used and their auxiliary functions can not be found:
!       Restrict to '-DABAQUS_CALLING' and check that this either sets 'ABQ_STD_CALLING' or 'ABQ_EXP_CALLING' internally.
! If an Abaqus simulation fails with "ERROR: NO VUMAT SUBROUTINE WAS SUPPLIED FOR THE USER MATERIAL NAMED ...", possible reasons are:
!    -> the correct prefix is missing in the material name.
!    -> the wrong interface is called. Check that ABAQUS_CALLING points to the correct macro (ABQ_STD_CALLING or ABQ_EXP_CALLING).
! If an Abaqus simulation fails with "Abaqus Error: Abaqus/Standard Analysis exited with an error" or
! "Abaqus Error: Abaqus/Explicit Packager exited with an error" without creating ".dat" or ".sta files:
!   -> One possibility is that neither "ABQ_STD_CALLING" nor "ABQ_EXP_CALLING" is defined and therefore
!      the xMat not propberly linked to the simulation.
! If a Fortran runtime error is thrown (e.g. with index of dimension above upper bound):
!    -> Possibly the calling program uses a different (order of) arguments when calling.
!       Check the current interface and compare it with the expected interface of the calling program.
!
! ---------------------------
! -   R e f e r e n c e s   -
! ---------------------------
! Based on the works of P.-A. von Wolffersdorff, A. Niemunis and I. Herle,
! inspired by the original Fortran scripts from A. Niemunis
! with modifications by M. Kelm (2002), K.-P. Mahutka (2004), T. Pucker (2011) and G. Qiu (2012)
!
! Main literature used on constitutive models:
!  - Bauer, E. (1996): "Calibration of a comprehensive hypoplastic model for granular materials". In: Soils and Foundations 36, p.13--26.
!  - Niemunis, A. (2003): "Extended hypoplastic models for soils", Habilitation. Ruhr-Universitaet Bochum.
!  - Schranz, F. (2018): "Applications and developments of Barodesy", PhD Thesis. University of Innsbruck.
!  - von Wolffersdorff, P.-A. (1996): "A hypoplastic relation for granular materials with a predefined limit state surface".
!       In: Mechanics of Cohesive-frictional Materials 1, p.215--271.
!
! Main literature used on solver and mathematics:
!  - Press, W. H., Teukolsky S. A., Vetterling W. T. and Flannery B. P. (2007): "Numerical Recipes". ISBN: 978-0-521-88068-6

! Main literature used on interfaces:
!  - SIMULIA Abaqus (2020): "User defined subroutines". https://help.3ds.com
!  - Mathworks Matlab (2019): "Components of a Fortran MEX File".
!       https://www.mathworks.com/help/matlab/matlab_external/components-of-fortran-mex-file.html
!  - Plaxis (2016): "User-defined soil models".
!       https://communities.bentley.com/products/geotech-analysis/w/plaxis-soilvision-wiki/45468/plaxis-user-defined-soil-models
!  - Tochnog (2021): "User's manual". https://www.tochnogprofessional.nl/manuals/page.php
!       and "README_UMAT_USER.txt" from https://www.tochnogprofessional.nl/download/tochnog_linux_64_bit_user_supplied.tar.gz
!
! Main literature used on Fortran:
!  - Metcalf, M., Reid. J., Cohen, M. (2011): "Modern Fortran explained". ISBN: 978-0-199-60141-7.
!
! The official reference for this routine is the PhD thesis
!  - Zobel, D. (2023): "Zur Implementierung hoeherwertiger Stoffmodelle am Beispiel der Hypoplastizitaet", Dissertation. TU Hamburg.
!
! Code is best viewed when compiled in LaTeX with listings package (140 chars wide), text in backticks interpreted as math and some visual
! adjustments (see code_printer.tex in the repository)


! Using preprocessor directives to fit routine to calling program. Abaqus and Matlab have to explicitly set one
! (here ABAQUS_CALLING and MATLAB_CALLING respectively). Matlab requires inclusion of a header file
#ifdef MATLAB_CALLING
#include "fintrf.h"
#endif
! Abaqus (2020) automatically introduces a preprocessor variable ABQ_FORTRAN, which might be used here instead of ABAQUS_CALLING
#ifdef ABAQUS_CALLING
! NOTE: Currently either ABQ_STD_CALLING has to be defined (when using Abaqus/Standard) or
!       ABQ_EXP_CALLING (when using Abaqus/Explicit or Abaqus/CEL).
#define ABQ_STD_CALLING
#endif
! Also make sure that CINTER is active when OCTAVE_CALLING is defined
#ifdef OCTAVE_CALLING
#define CINTER
#endif


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
   !                                                                 ! increment
   logical, parameter :: setting_restrict_initial_substep = .True.   ! Uses setting_initial_substep_scale for the first step (or all steps
   !                                                                 ! if integration method is non-adaptive) if set to .True. Otherwise
   !                                                                 ! the given dt of the calling program is used (resulting in a single
   !                                                                 ! step for non-adaptive methods)
   real(dp), parameter :: setting_initial_substep_scale = 0.00001_dp ! If enforced with setting_restrict_initial_substep the minimum of
   !                                                                 ! this and the given dt will be used for the (first) time increment
   real(dp), parameter :: setting_stepsize_scaling_safety = 0.9      ! Safety factor applied to the estimation of new time increments
   !                                                                 ! during the integration loop (see also Press, 1997, p. 712)
   real(dp), parameter :: setting_max_rel_error = 0.001_dp           ! Maximum relative error for adaptive integration methods. The error
   !                                                                 ! is scaled componentwise first and then the maximum absolute error
   !                                                                 ! has to be smaller than this value (see also Integrate subroutine in
   !                                                                 ! the solver module)
   integer, parameter :: setting_max_integration_refinements = 20    ! When the error is too large, an integration step is repeated with a
   !                                                                 ! smaller time increment (and checked again). Stop if those checks
   !                                                                 ! fail more often than this value in a row without any accepted step

   ! --- Internal memory reservation for Matlab (mexFunction) and Plaxis (USER_MOD). Might need to be adjusted for new constitutive models
   integer, parameter :: setting_num_statevariables = 20             ! Maximum number of state variables for each constitutive model
   integer, parameter :: setting_max_mat_params = 16                 ! Maximum number of provided material parameters
   integer, parameter :: setting_max_internal_states = &             ! Internal states include stresses (+jac) and statevariables (+jac)
                                                       (6+1)*6 + (6+1)*setting_num_statevariables

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


   ! --------------------------------------------------------------- !
   pure subroutine Estimate_Components(nel, num_dimensions, num_shear)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel                                     ! The number of elements nel has to be between 3 and 6
      integer, intent(out) :: num_dimensions, num_shear
      ! ------------------------------------------------------------ !
      num_dimensions = 0
      num_shear = 0

      if (nel == 3) then                                             ! Either two direct and one shear component
         num_dimensions = 2
         num_shear = 1
      else if ((nel > 3) .and. (nel <= 6)) then                      ! Or three direct and up to three shear components for 3 < nel <= 6
         num_dimensions = 3
         num_shear = nel - num_dimensions
      end if
   end subroutine Estimate_Components


   ! --------------------------------------------------------------- !
   function Check_Input_Dimensions(num_dimensions, num_shear)
   ! --------------------------------------------------------------- ! Only support 2D and 3D symmetric input
      integer, intent(in) :: num_dimensions, num_shear
      logical :: Check_Input_Dimensions
      ! ------------------------------------------------------------ !
      Check_Input_Dimensions = .False.

      if (num_dimensions == 2) then
         if (num_shear == 1) then
            Check_Input_Dimensions = .True.
         end if
      else if (num_dimensions == 3) then
         if ((num_shear >= 1) .and. (num_shear <= 3)) then
            Check_Input_Dimensions = .True.
         end if
      end if

      if (Check_Input_Dimensions) then
         global_num_direct_components = num_dimensions               ! Changing global variables makes this function impure
         global_num_shear_components = num_shear
      end if
   end function Check_Input_Dimensions


   ! --------------------------------------------------------------- !
   elemental function Is_Nan(number)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: number
      logical :: Is_Nan
      ! ------------------------------------------------------------ !
      if (number /= number) then                                     ! Only true for NaN
         Is_Nan = .True.
      else
         Is_Nan = .False.
      end if
   end function Is_Nan


   ! --------------------------------------------------------------- !
   pure function C_To_Fortran_String(length, c_string) result(f_string)
   ! --------------------------------------------------------------- !
      use iso_c_binding, only: c_char, c_null_char

      integer, intent(in) :: length
      character(kind=c_char, len=1), dimension(length), intent(in) :: c_string
      character(len=length) :: f_string
      ! ------------------------------------------------------------ !
      integer :: idx

      f_string = " "
      do idx = 1, length
         if (c_string(idx) == c_null_char) then
            exit
         else
            f_string(idx:idx) = c_string(idx)
         end if
      end do
   end function C_To_Fortran_String


   ! --------------------------------------------------------------- !
   pure function Uppercase(text) result(upper)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: text
      character(len=len(text)) :: upper
      ! ------------------------------------------------------------ !
      character(len=26), parameter :: low_alph = 'abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: up_alph  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer :: idx, pos

      upper = text
      do idx = 1, len(text)
         pos = index(low_alph, text(idx:idx))
         if (pos /= 0) then
            upper(idx:idx) = up_alph(pos:pos)
         end if
      end do
   end function Uppercase


   ! --------------------------------------------------------------- !
   subroutine Write_Message(message)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: message
      ! ------------------------------------------------------------ !
      character(len=len(identifier_log_start) + len(message)) :: InfoMessage

      InfoMessage = identifier_log_start // message

      ! For Abaqus calls either stdb_abqerr (Abaqus/Standard) or xplb_abqerr (Abaqus/Explicit) is
      ! used. Besides a printf-ed message string, an int (array), a real (array) and a char (array)
      ! are expected - See Abaqus user routine guide
#ifdef ABQ_STD_CALLING
      call stdb_abqerr(1, InfoMessage, 0, 0.0, ' ')                  ! Output on Abaqus message log channel (1)
#else
#ifdef ABQ_EXP_CALLING
      call xplb_abqerr(1, InfoMessage, 0, 0.0, ' ')                  ! Output on Abaqus message log channel (1)
#else
#ifdef MATLAB_CALLING
      call mexPrintf(InfoMessage)
#else
#ifdef PLAXIS_DLL
      write(1, *) InfoMessage                                        ! Output on Plaxis log channel (1)
#else
      write(0, *) InfoMessage                                        ! Output on default log channel (0)
#endif
#endif
#endif
#endif
   end subroutine Write_Message


   ! --------------------------------------------------------------- !
   subroutine Write_Warning(message)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: message
      ! ------------------------------------------------------------ !
      character(len=len(identifier_log_start) + len(message)) :: WarningMessage

      WarningMessage = identifier_log_start // message

      ! For Abaqus calls either stdb_abqerr (Abaqus/Standard) or xplb_abqerr (Abaqus/Explicit) is
      ! used. Besides a printf-ed message string, an int (array), a real (array) and a char (array)
      ! are expected - See Abaqus user routine guide
#ifdef ABQ_STD_CALLING
      call stdb_abqerr(-1, WarningMessage, 0, 0.0, ' ')              ! Output on Abaqus warning log channel (-1)
#else
#ifdef ABQ_EXP_CALLING
      call xplb_abqerr(-1, WarningMessage, 0, 0.0, ' ')              ! Output on Abaqus warning log channel (-1)
#else
#ifdef MATLAB_CALLING
      call mexWarnMsgTxt(WarningMessage)
#else
#ifdef PLAXIS_DLL
      write(1, *) WarningMessage                                     ! Output on Plaxis log channel (1)
#else
      write(0, *) WarningMessage                                     ! Output on default log channel (0)
#endif
#endif
#endif
#endif
   end subroutine Write_Warning


   ! --------------------------------------------------------------- !
   subroutine Write_Error_And_Exit(message)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: message
      ! ------------------------------------------------------------ !
      character(len=len(identifier_log_start) + len(message)) :: ErrorMessage

      ErrorMessage = identifier_log_start // message

      ! For Abaqus calls either stdb_abqerr (Abaqus/Standard) or xplb_abqerr (Abaqus/Explicit) is
      ! used. Besides a printf-ed message string, an int (array), a real (array) and a char (array)
      ! are expected - See Abaqus user routine guide
#ifdef ABQ_STD_CALLING
      call stdb_abqerr(-2, ErrorMessage, 0, 0.0, ' ')                ! Output on Abaqus error log channel (-2)
      call xit()
#else
#ifdef ABQ_EXP_CALLING
      call xplb_abqerr(-2, ErrorMessage, 0, 0.0, ' ')                ! Output on Abaqus error log channel (-2)
      call xplb_exit()
#else
#ifdef MATLAB_CALLING
      call mexErrMsgTxt(ErrorMessage)
#else
#ifdef PLAXIS_DLL
      write(1, *) ErrorMessage                                       ! Output on Plaxis log channel (1)
      stop 'Execution aborted'
#else
#ifdef OCTAVE_CALLING
      write(0, *) ErrorMessage
      call xStopx('Execution aborted')
#else
      write(0, *) ErrorMessage                                       ! Output on default log channel (0)
      stop 'Execution aborted'
#endif
#endif
#endif
#endif
#endif
   end subroutine Write_Error_And_Exit
end module General_Settings


! ==================================================================================================================== !
module Debug
   use General_Settings, only: dp
   implicit none

   integer, parameter :: onenumberlength = 12
   integer, parameter :: mat_identifierlength = 6
   integer, parameter :: tens_identifierlength = 12
   character(len=6), parameter :: submatrixname = '     M'
   character(len=6), parameter :: debug_format = 'f12.5 '            ! character length of debug_format and
   character(len=6), parameter :: alternative_format = 'es12.3'      ! alternative_format must be the same

   interface Formatval
      module procedure Format_logical, Format_int_dim0, Format_dp_dim0, Format_dp_dim1, Format_dp_dim2, Format_dp_dim4
   end interface

   private
   public Formatval


   contains


   ! --------------------------------------------------------------- !
   elemental function Number_To_String(number)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan
      !
      real(dp), intent(in) :: number
      character(len=onenumberlength+2) :: Number_To_String
      ! ------------------------------------------------------------ !
      character(len=len(debug_format)+2) :: formatstring

      Number_To_String(1:onenumberlength+2) = repeat(' ', onenumberlength+2)
      Number_To_String(onenumberlength+1:onenumberlength+1) = ','

      if (Is_Nan(number)) then
         write(Number_To_String(onenumberlength-2:onenumberlength), '(a)') 'NaN'
      else                                                           ! Use fixed format for numbers "close" to zero
         if ((abs(number) > 1000.0_dp) .or. (abs(number) < 0.001_dp)) then
            if (abs(number) < setting_epsilon) then
               write(Number_To_String(onenumberlength-2:onenumberlength), '(a)') '0.0'
            else
               write(formatstring, '(a, a, a)') '(', alternative_format, ')'
               write(Number_To_String(1:onenumberlength), fmt=formatstring) number
            end if
         else
            write(formatstring, '(a, a, a)') '(', debug_format, ')'
            write(Number_To_String(1:onenumberlength), fmt=formatstring) number
         end if
      end if
   end function Number_To_String


   ! --------------------------------------------------------------- !
   pure function Format_logical(name, bool)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      logical, intent(in) :: bool
      character(len=len(name)+6) :: Format_logical
      ! ------------------------------------------------------------ !
      if (bool) then
         write(Format_logical, '(a, a)') name, ':  Yes'
      else
         write(Format_logical, '(a, a)') name, ':  No '
      end if
   end function Format_logical


   ! --------------------------------------------------------------- !
   pure function Format_int_dim0(name, number)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      integer, intent(in) :: number
      character(len=len(name)+2 + 6) :: Format_int_dim0
      ! ------------------------------------------------------------ !
      write(Format_int_dim0, '(a, a, i6)') name, '  ', number
   end function Format_int_dim0


   ! --------------------------------------------------------------- !
   pure function Format_dp_dim0(name, number)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: number
      character(len=len(name)+2 + onenumberlength+2) :: Format_dp_dim0
      ! ------------------------------------------------------------ !
      write(Format_dp_dim0, '(a, a, a)') name, '  ', Number_To_String(number)
      Format_dp_dim0(len(Format_dp_dim0)-1:len(Format_dp_dim0)-1) = ' '
   end function Format_dp_dim0


   ! --------------------------------------------------------------- !
   pure function Format_dp_dim1(name, vector)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:), intent(in) :: vector
      character(len=len(name)+2 + (onenumberlength+2)*size(vector)) :: Format_dp_dim1
      ! ------------------------------------------------------------ !
      character(len=11) :: formatstring

      write(formatstring, '(a, i2, a)') '(a, a, ', size(vector), 'a)'
      write(Format_dp_dim1, formatstring) name, '  ', Number_To_String(vector)
      Format_dp_dim1(len(Format_dp_dim1)-1:len(Format_dp_dim1)-1) = ' '
      Format_dp_dim1(len(Format_dp_dim1):len(Format_dp_dim1)) = char(10)
   end function Format_dp_dim1


   ! --------------------------------------------------------------- !
   pure function Format_dp_dim2(name, matrix)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:, :), intent(in) :: matrix
      character(len=(len(name)+2 + mat_identifierlength+2 + &
         (onenumberlength+2)*size(matrix, 1))*size(matrix, 2)) :: Format_dp_dim2
      ! ------------------------------------------------------------ !
      character(len=mat_identifierlength) :: identifier
      integer :: zeilenlaenge, idx

      zeilenlaenge = len(name)+2 + mat_identifierlength+2 + (onenumberlength+2)*size(matrix, 1)
      identifier = '(:, :)'
      do idx = 1, size(matrix, 2)
         write(identifier(2:2), '(i1)') idx
         write(Format_dp_dim2((idx-1)*zeilenlaenge+1:idx*zeilenlaenge), '(a, a, a, a)') &
            name, identifier, '  ', Format_dp_dim1(name='', vector=matrix(idx, :))
      end do
   end function Format_dp_dim2


   ! --------------------------------------------------------------- !
   pure function Format_dp_dim4(name, tensor)
   ! --------------------------------------------------------------- !
      character(len=*), intent(in) :: name
      real(dp), dimension(:, :, :, :), intent(in) :: tensor
      character(len=((len(submatrixname)+2+mat_identifierlength+2 + (onenumberlength+2)*size(tensor, 4)) &
         *size(tensor, 3) + len(name)+tens_identifierlength+1)*size(tensor, 2)*size(tensor, 1)) :: Format_dp_dim4
      ! ------------------------------------------------------------ !
      character(len=tens_identifierlength) :: identifier
      integer :: matrixlaenge, idx1, idx2, position

      matrixlaenge = (len(submatrixname)+2+mat_identifierlength+2 + (onenumberlength+2)*size(tensor, 4)) &
                   * size(tensor, 3) + len(name)+tens_identifierlength+1
      identifier = '(:, :, :, :)'
      do idx1 = 1, size(tensor, 1)
         write(identifier(2:2), '(i1)') idx1
         do idx2 = 1, size(tensor, 2)
            position = ((idx1-1)*size(tensor, 2) + (idx2-1))*matrixlaenge
            write(identifier(5:5), '(i1)') idx2

            write(Format_dp_dim4(position+1:position+matrixlaenge), '(a, a, a, a)') &
               name, identifier, char(10), Format_dp_dim2(name=submatrixname, matrix=tensor(idx1, idx2, :, :))
         end do
      end do
   end function Format_dp_dim4
end module Debug


! ==================================================================================================================== !
module Math_Operations
   use General_Settings, only: dp, setting_epsilon, global_num_direct_components, global_num_shear_components
   implicit none

   real(dp), parameter :: const_root2 = sqrt(2.0_dp)
   real(dp), parameter :: const_root3 = sqrt(3.0_dp)
   real(dp), parameter :: const_root6 = sqrt(6.0_dp)
   !
   ! For 6x1 vectors and 6x6 matrices the tensorial index pairs `(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 3)` in both dimensions are used
   real(dp), parameter, dimension(6) :: const_identity2d = &         ! 2dim identity matrix `\mathbf{I}` where `I_{ij} = \delta_{ij}`
      [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]               ! with Kronecker-symbol `\delta_{ij} = 1` for `i=j` and 0 else
   real(dp), parameter, dimension(6, 6) :: const_identity4d_sym = &  ! 4dim symmetric identity tensor `\mathcal{I}^\mathrm{sym}` with `I^\mathrm{sym}_{ijkl} = \frac{1}{2}\left(\delta_{ik}\otimes\delta_{jl} + \delta_{il}\otimes\delta_{jk}\right)`
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &     ! which returns the symmetric part of the matrix: `\mathcal{I}^\mathrm{sym}:\mathbf{A} = \mathrm{sym}\left(\mathbf{A}\right)`
               0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [6, 6])
   real(dp), parameter, dimension(6, 6) :: const_identity4d_tr = &   ! 4dim identity tensor representing `\mathcal{I}^\mathrm{tr}` with `I^\mathrm{tr}_{ijkl} = \delta_{ij}\otimes\delta_{kl}`
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &     ! which maps the trace of a matrix to each element of a
               1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &     ! diagonal matrix:  `\mathcal{I}^\mathrm{tr}:\mathbf{A} = \tr{(\mathbf{A})}\cdot\delta_{ii}`
               1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [6, 6])
   integer, parameter, dimension(3, 3) :: ref_elements = reshape([ &
      1, 4, 6, &
      4, 2, 5, &
      6, 5, 3], [3, 3])
   integer, parameter, dimension(2, 6) :: ref_indices = reshape([ &
      1, 1,   2, 2,   3, 3,   1, 2,   2, 3,   1, 3], [2, 6])

   private
   public const_root2, const_root3, const_root6, &
          const_identity2d, const_identity4d_sym, const_identity4d_tr, &
          Nonzero_Division, Trace, Dimensionless, Deviatoric_Part, Norm, Determinant, Squared, &
          Inverse_Tensor, Double_Contraction22, Double_Contraction42, Double_Contraction44, Dyadic_Product22, &
          Tensor_Partialtrace, Matrix_Rotation, Matrix_Exponential, Sqrt_Of_Sum_Of_Squares, &
          Is_Diagonal, Value_In_Interval, Abort_If_Not_In_Interval, Set_Element_In_Tensor, &
          Get_Matrix_From_Tensorpart, Set_Matrix_In_Tensorpart, Get_Elements_From_Matrixlist, &
          Set_Elements_In_Matrixlist, Ref_Index, Perturbate, Vec9_To_Mat, Mat_To_Vec9, &
          Pack_States, Unpack_States, Import_Matrix, Export_Matrix, Export_Tensor


   contains


   ! --------------------------------------------------------------- !
   elemental function Nonzero_Division(val, fac)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: val
      real(dp), intent(in) :: fac
      real(dp) :: Nonzero_Division
      ! ------------------------------------------------------------ !
      if (abs(fac) < setting_epsilon) then
         Nonzero_Division = val
      else
         Nonzero_Division = val/fac
      end if
   end function Nonzero_Division


   ! --------------------------------------------------------------- !
   pure function Trace(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: Trace
      ! ------------------------------------------------------------ !
      Trace = sum(vec6(1:global_num_direct_components))              ! Trace `\tr{(\mathbf{M})} = \sum_{i=1}^n M_{ii}`
   end function Trace


   ! --------------------------------------------------------------- !
   pure function Dimensionless(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: Dimensionless
      ! ------------------------------------------------------------ !
      Dimensionless = Nonzero_Division(val=vec6, fac=Trace(vec6))    ! Dimensionless `\hat{\mathbf{M}} = \frac{\mathbf{M}}{\tr{(\mathbf{M})}}`
   end function Dimensionless


   ! --------------------------------------------------------------- !
   pure function Deviatoric_Part(vec6) result(dev6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: dev6
      ! ------------------------------------------------------------ !
      dev6 = vec6                                                    ! Deviatoric part (3D) `\mathbf{M}^\mathrm{dev} = \mathbf{M} - \frac{\tr{(\mathbf{M})}}{n}\mathbf{I}`
      associate(ndi => global_num_direct_components)
         dev6(1:ndi) = dev6(1:ndi) - Trace(vec6)/ndi
      end associate
   end function Deviatoric_Part


   ! --------------------------------------------------------------- !
   pure function Norm(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: Norm
      ! ------------------------------------------------------------ !
      Norm = sqrt(Double_Contraction22(vec6, vec6))                  ! Matrix 2-Norm (Frobenius norm) `||\mathbf{M}|| = \sqrt{\sum_i\sum_j|M_{ij}|^2}`
   end function Norm


   ! --------------------------------------------------------------- !
   pure function Determinant(vec6) result(det6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp) :: det6
      ! ------------------------------------------------------------ !
      if (global_num_direct_components == 3) then
         det6 = vec6(1)*(vec6(2)*vec6(3) - vec6(5)*vec6(5)) &
              + vec6(4)*(vec6(5)*vec6(6) - vec6(3)*vec6(4)) &
              + vec6(6)*(vec6(4)*vec6(5) - vec6(2)*vec6(6))
      else
         det6 = vec6(1)*vec6(2) - vec6(4)*vec6(4)
      end if
   end function Determinant


   ! --------------------------------------------------------------- !
   pure function Squared(vec6) result(sq_vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: sq_vec6
      ! ------------------------------------------------------------ !
      sq_vec6(1) = vec6(1)**2 + vec6(4)**2 + vec6(6)**2              ! Square `\mathbf{M}^2 = \mathbf{M}\cdot \mathbf{M}`
      sq_vec6(2) = vec6(4)**2 + vec6(2)**2 + vec6(5)**2              ! (square of a symmetric matrix is also symmetric)
      sq_vec6(3) = vec6(6)**2 + vec6(5)**2 + vec6(3)**2
      sq_vec6(4) = vec6(4)*(vec6(1) + vec6(2)) + vec6(5)*vec6(6)
      sq_vec6(5) = vec6(5)*(vec6(2) + vec6(3)) + vec6(4)*vec6(6)
      sq_vec6(6) = vec6(6)*(vec6(1) + vec6(3)) + vec6(4)*vec6(5)
   end function Squared


   ! --------------------------------------------------------------- !
   pure subroutine LU_Decomposition(matrix, nel, trans_mat, success)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! Determine the LU-decomposition of matrix `\mathbf{A}`, so that `\mathbf{P}\mathbf{A} = \mathbf{L}\mathbf{U}`
      real(dp), dimension(nel, nel), intent(out) :: trans_mat        ! Save `\mathbf{L}` and `\mathbf{U}` in matrix and the permutation `\mathbf{P}` in trans_mat
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), parameter, dimension(2, 2) :: permutation = reshape([ &
         0.0_dp, 1.0_dp, &
         1.0_dp, 0.0_dp], [2, 2])
      integer :: idx, jdx, idx_max

      success = .True.
      trans_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])

      do idx = 1, nel-1
         ! Find the largest element for (partial) pivoting and adjust transformation matrix
         idx_max = idx + maxloc(abs(matrix(idx:nel, idx)), dim=1) - 1
         if (idx_max /= idx) then
            trans_mat([idx, idx_max], :) = matmul(permutation, trans_mat([idx, idx_max], :))
            matrix([idx, idx_max], :) = matrix([idx_max, idx], :)
         end if

         if (abs(matrix(idx, idx)) < setting_epsilon) then
            ! If value on diagonal (maximum value of matrix(idx:nel, idx)) is zero, the matrix is singular
            success = .False.
            cycle
         end if

         matrix(idx+1:nel, idx) = matrix(idx+1:nel, idx)/matrix(idx, idx)
         matrix(idx+1:nel, idx+1:nel) = matrix(idx+1:nel, idx+1:nel) &
                                      - matmul(reshape(matrix(idx+1:nel, idx), [nel-idx, 1]), &
                                               reshape(matrix(idx, idx+1:nel), [1, nel-idx]))
      end do
   end subroutine LU_Decomposition


   ! --------------------------------------------------------------- !
   pure subroutine LU_Solve(lu_mat, nel, trans_mat, b_vec)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: lu_mat            ! Use LU decomposition of matrix `\mathbf{A}` saved in lu_mat and trans_mat
      real(dp), dimension(nel, nel), intent(in) :: trans_mat         ! representing `\mathbf{P}\mathbf{A} = \mathbf{L}\mathbf{U}` to calculate `\mathbf{A}x = b` with `\mathbf{L}y = \mathbf{P}b` and `\mathbf{U}x = y`.
      real(dp), dimension(nel), intent(inout) :: b_vec
      ! ------------------------------------------------------------ !
      integer :: idx

      ! All results are directly be written to b_mat, but association is used for code clarity
      associate(pb_vec => b_vec, y_vec => b_vec, x_vec => b_vec)
         pb_vec = matmul(trans_mat, pb_vec)                          ! `b^\prime = \mathbf{P}b`

         do idx = 2, nel                                             ! Solve for `y` in `\mathbf{L}y = \mathbf{P}b` by forward substitution
            y_vec(idx) = pb_vec(idx) &                               ! `y_i = \left(b_i^\prime \sum_{j=1}^{i-1} L_{ij}y_j\right) \frac{1}{L_{ii}}` with `L_{ii} = 1` and therefore `y_1 = b_1^\prime`
                       - dot_product(lu_mat(idx, 1:idx-1), y_vec(1:idx-1))
         end do

         x_vec(nel) = y_vec(nel)/lu_mat(nel, nel)                    ! Solve for `x` in `\mathbf{U}x = y` by back substitution
         do idx = nel-1, 1, -1                                       ! `x_i = \left(y_i - \sum_{j=i+1}^{n} U_{ij} x_j\right) \frac{1}{U_{ii}}`
            x_vec(idx) = (y_vec(idx) - dot_product(lu_mat(idx, idx+1:nel), x_vec(idx+1:nel)))/lu_mat(idx, idx)
         end do
         ! By association, all values determined for x_vec are already stored in b_vec
      end associate
   end subroutine LU_Solve


   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Internal(matrix, inv_matrix, nel, success)
   ! --------------------------------------------------------------- ! Calculates the inverse of a nel by nel matrix if it is not singular
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: matrix
      real(dp), dimension(nel, nel), intent(out) :: inv_matrix
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: lu_mat, p_mat
      integer :: idx, jdx

      lu_mat = matrix
      inv_matrix = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])

      ! Calculates LU-decomposition of `\mathbf{A}` and solves `Ax_j = \delta_j` for each column vector `\delta_j` of the identity matrix to get `\left[x_1| \dots |x_n\right] = \mathbf{A}^{-1}`
      call LU_Decomposition(matrix=lu_mat, nel=size(matrix, 1), trans_mat=p_mat, success=success)
      if (success) then
         do idx = 1, nel
            call LU_Solve(lu_mat=lu_mat, nel=nel, trans_mat=p_mat, b_vec=inv_matrix(:, idx))
         end do
      end if
   end subroutine Inverse_Internal


   ! --------------------------------------------------------------- !
   pure subroutine Inverse_Tensor(tensor, inv_tensor, success)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: tensor
      real(dp), dimension(6, 6), intent(out) :: inv_tensor
      logical, intent(out) :: success
      ! ------------------------------------------------------------ !
      real(dp), dimension(global_num_direct_components+global_num_shear_components, &
         global_num_direct_components+global_num_shear_components) :: comp_matrix, inv_comp_matrix

      ! To successfully find an inverse, the matrix has to have full rank. Since all calculations are done with vec6 or mat66,
      ! for inversion all rows and columns added for simplicity have to be temporarily removed.
      ! If the inversion fails, success is set to .False. and the calling function has to deal with it.
      associate(ndi => global_num_direct_components, nshr => global_num_shear_components)
         comp_matrix(1:ndi, 1:ndi) = tensor(1:ndi, 1:ndi)
         comp_matrix((ndi + 1):(ndi + nshr), 1:ndi) = tensor(4:(3 + nshr), 1:ndi)
         comp_matrix(1:ndi, (ndi + 1):(ndi + nshr)) = tensor(1:ndi, 4:(3 + nshr))
         comp_matrix((ndi + 1):(ndi + nshr), (ndi + 1):(ndi + nshr)) = tensor(4:(3 + nshr), 4:(3 + nshr))

         call Inverse_Internal(matrix=comp_matrix, inv_matrix=inv_comp_matrix, nel=ndi+nshr, success=success)

         inv_tensor = 0.0_dp
         inv_tensor(1:ndi, 1:ndi) = inv_comp_matrix(1:ndi, 1:ndi)
         inv_tensor(4:(3 + nshr), 1:ndi) = inv_comp_matrix((ndi + 1):(ndi + nshr), 1:ndi)
         inv_tensor(1:ndi, 4:(3 + nshr)) = inv_comp_matrix(1:ndi, (ndi + 1):(ndi + nshr))
         inv_tensor(4:(3 + nshr), 4:(3 + nshr)) = inv_comp_matrix((ndi + 1):(ndi + nshr), (ndi + 1):(ndi + nshr))
      end associate
   end subroutine Inverse_Tensor


   ! --------------------------------------------------------------- !
   pure function Double_Contraction22(vec6a, vec6b) result(vec_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6a, vec6b
      real(dp) :: vec_out
      ! ------------------------------------------------------------ !
      vec_out = dot_product(vec6a, [vec6b(1:3), 2.0_dp*vec6b(4:6)])  ! Double Contraction `\mathbf{M}_{A} : \mathbf{M}_{B} = \sum_{ij} M_{A,ij}\cdot M_{B,ij}`
   end function Double_Contraction22


   ! --------------------------------------------------------------- !
   pure function Double_Contraction42(mat66, vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(6) :: Double_Contraction42
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: temp_vec

      temp_vec(1:3) = vec6(1:3)
      temp_vec(4:6) = 2.0_dp*vec6(4:6)
      Double_Contraction42 = matmul(mat66, temp_vec)                 ! Double contraction `\mathcal{T} : \mathbf{M} = \sum_{kl} T_{ijkl}\cdot M_{kl}`
   end function Double_Contraction42


   ! --------------------------------------------------------------- !
   pure function Double_Contraction44(mat66a, mat66b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66a, mat66b
      real(dp), dimension(6, 6) :: Double_Contraction44
      ! ------------------------------------------------------------ !
      real(dp), dimension(6, 6) :: temp_mat66a

      temp_mat66a(:, 1:3) = mat66a(:, 1:3)
      temp_mat66a(:, 4:6) = 2.0_dp*mat66a(:, 4:6)
      Double_Contraction44 = matmul(temp_mat66a, mat66b)             ! Double contraction `\mathcal{T} : \mathcal{F} = \sum_{kl} T_{ijkl}\cdot F_{klmn}`
   end function Double_Contraction44


   ! --------------------------------------------------------------- !
   pure function Dyadic_Product22(vec6a, vec6b)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6a, vec6b
      real(dp), dimension(6, 6) :: Dyadic_Product22
      ! ------------------------------------------------------------ !
      Dyadic_Product22 = matmul(reshape(vec6a, [6, 1]), &            ! Dyadic product `T_{ijkl} = M_{A,ij}\cdot M_{B,kl}`
                                reshape(vec6b, [1, 6]))
   end function Dyadic_Product22


   ! --------------------------------------------------------------- !
   pure function Tensor_Partialtrace(mat66)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: mat66
      real(dp) :: Tensor_Partialtrace
      ! ------------------------------------------------------------ !
      Tensor_Partialtrace = mat66(1, 1) + mat66(2, 2) + mat66(3, 3)
   end function Tensor_Partialtrace


   ! --------------------------------------------------------------- !
   pure function Matrix_Rotation(vec6, rot_mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(3, 3), intent(in) :: rot_mat33
      real(dp), dimension(6) :: Matrix_Rotation
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: M_mat

      M_mat = Vec_To_Mat33(vec6=vec6)
      M_mat = matmul(rot_mat33, matmul(M_mat, transpose(rot_mat33))) ! Rotation `\mathbf{D}\cdot \mathbf{X}\cdot \mathbf{D}^T` of symmetric matrix `\mathbf{X}`
      Matrix_Rotation = Mat33_To_Vec(mat33=M_mat)                    ! results in another symmetric matrix
   end function Matrix_Rotation


   ! --------------------------------------------------------------- !
   elemental function Sqrt_Of_Sum_Of_Squares(a, b)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: a, b
      real(dp) :: Sqrt_Of_Sum_Of_Squares
      ! ------------------------------------------------------------ !
      real(dp) :: abs_value_a, abs_value_b

      abs_value_a = abs(a)
      abs_value_b = abs(b)

      if (abs_value_a > abs_value_b) then
         Sqrt_Of_Sum_Of_Squares = abs_value_a*sqrt(1.0_dp + (abs_value_b/abs_value_a)**2)
      else if (abs_value_b < setting_epsilon) then
         Sqrt_Of_Sum_Of_Squares = 0.0_dp
      else
         Sqrt_Of_Sum_Of_Squares = abs_value_b*sqrt(1.0_dp + (abs_value_a/abs_value_b)**2)
      end if
   end function Sqrt_Of_Sum_Of_Squares


   ! --------------------------------------------------------------- !
   pure subroutine Tridiagonal_Matrix(matrix, nel, trans_mat)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! Computes a tridiagonal matrix and an assembly of transformation
      real(dp), dimension(nel, nel), intent(out) :: trans_mat        ! matrices as described in section 11.2 of Press et al (1997)
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: ident_mat, householder_mat    ! The explicit creation of the householder matrices (as done here)
      real(dp), dimension(nel) :: u_vec                              ! is only justifiable for small matrices (nel << 10)
      real(dp) :: H_mat, sigma, scaling_factor
      integer :: idx, jdx, ldx

      ident_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])
      trans_mat = ident_mat                                          ! `\mathbf{Q}_0 = \mathbf{I}`

      do idx = nel, 3, -1
         ldx = idx - 1
         scaling_factor = sum(abs(matrix(idx, 1:ldx)))

         if (scaling_factor > 0.0_dp) then
            sigma = dot_product(matrix(idx, 1:ldx)/scaling_factor, matrix(idx, 1:ldx)/scaling_factor)
            sigma = sign(sqrt(sigma), matrix(idx, ldx))

            u_vec = 0.0_dp
            u_vec(1:ldx) = matrix(idx, 1:ldx)/scaling_factor
            u_vec(ldx) = u_vec(ldx) + sigma
            H_mat = 0.5_dp*dot_product(u_vec(1:ldx), u_vec(1:ldx))   ! (11.2.4) of Press et al (1997): `H = \frac{1}{2}|u|^2`

            householder_mat = ident_mat                              ! (11.2.3) of Press et al (1997): `\mathbf{P}_{n+1} = \mathbf{I} - \frac{u u^T}{H}`
            householder_mat(1:ldx, 1:ldx) = householder_mat(1:ldx, 1:ldx) &
                                          - matmul(reshape(u_vec(1:ldx)/H_mat, [ldx, 1]), &
                                                   reshape(u_vec(1:ldx), [1, ldx]))
            matrix = matmul(householder_mat, matmul(matrix, transpose(householder_mat)))
            trans_mat = matmul(trans_mat, householder_mat)           ! (11.2.17) of Press et al (1997): `\mathbf{Q}_{n+1} = \mathbf{Q}_n \mathbf{P}_{n+1}`
         end if
      end do

      ! Set all elements above/beneath subdiagonal to exactly zero
      do idx = 1, nel
         do jdx = 1, nel
            if (abs(jdx-idx) > 1) then
               matrix(jdx, idx) = 0.0_dp
            end if
         end do
      end do
   end subroutine Tridiagonal_Matrix


   ! --------------------------------------------------------------- !
   pure function Givens_CS(refval, zeroval)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: refval, zeroval
      real(dp), dimension(2, 2) :: Givens_CS
      ! ------------------------------------------------------------ !
      real(dp) :: tau, cos_theta, sin_theta

      if (abs(zeroval) < setting_epsilon) then
         cos_theta = 1.0_dp
         sin_theta = 0.0_dp
      else if (abs(zeroval) > abs(refval)) then
         tau = -refval/zeroval;
         sin_theta = 1.0_dp/sqrt(1.0_dp + tau**2)
         cos_theta = tau*sin_theta
      else
         tau = -zeroval/refval;
         cos_theta = 1.0_dp/sqrt(1.0_dp + tau**2)
         sin_theta = tau*cos_theta
      end if

      Givens_CS = reshape([cos_theta, sin_theta, -sin_theta, cos_theta], [2, 2])
   end function Givens_CS


   ! --------------------------------------------------------------- !
   pure subroutine QL_Decomposition(matrix, nel, eigenvec_mat, success)
   ! --------------------------------------------------------------- ! Use this function for symmetric tridiagonal matrices only.
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(inout) :: matrix         ! QL decomposition with implicit shifts of matrix as described in
      real(dp), dimension(nel, nel), intent(out) :: eigenvec_mat     ! section 11.3 of Press et al (1997) returning the eigenvalues in
      logical, intent(out) :: success                                ! matrix (not sorted), the eigenvectors and a success flag
      ! ------------------------------------------------------------ !
      integer, parameter :: max_iterations = 30
      integer :: idx, jdx, kdx, idx_start, idx_end, iter
      real(dp) :: diff, mu, refval, zeroval, temp_sum
      real(dp), dimension(2, 2) :: givens
      real(dp), parameter, dimension(2, 2) :: identity_2d = reshape([ &
         1.0_dp, 0.0_dp, &
         0.0_dp, 1.0_dp], [2, 2])
      logical :: extract_block

      eigenvec_mat = reshape([(1.0_dp, (0.0_dp, idx = 1, nel), jdx = 1, nel-1), 1.0_dp], [nel, nel])
      success = .False.

      do idx_start = 1, nel
         iter = 0
         do
            iter = iter+1
            if (iter == max_iterations) then
               ! Too many iterations: Cancel the decomposition procedure
               return
            end if

            extract_block = .False.
            ! Look for a single small subdiagonal element to split the matrix in two blocks.
            ! Focus on the top left block starting from (idx_start, idx_start)
            do idx_end = idx_start, nel-1
               temp_sum = abs(matrix(idx_end, idx_end)) + abs(matrix(idx_end+1, idx_end+1))
               ! If adding the subdiagonal element to temp_sum is smaller than the next representable
               ! value in the used precision, then it is negligible and the matrix can be split here
               if (nearest(temp_sum, 1.0_dp) >= temp_sum + abs(matrix(idx_end+1, idx_end))) then
                  extract_block = .True.
                  exit
               end if
            end do

            if (.not. extract_block) then
               idx_end = nel
            end if

            if (idx_end == idx_start) then
               exit
            end if

            ! The first rotation is special while following (idx_end-1-idx_start) rotations "chase"
            ! the only non-zero off-subdiagonal element out of the matrix
            diff = (matrix(idx_start+1, idx_start+1) - matrix(idx_start, idx_start))/2.0_dp
            mu = matrix(idx_start, idx_start) - matrix(idx_start+1, idx_start)**2 &
               / (diff + sign(Sqrt_Of_Sum_Of_Squares(diff, matrix(idx_start+1, idx_start)), diff))
            refval = matrix(idx_end, idx_end) - mu
            zeroval = matrix(idx_end, idx_end-1)

            do kdx = idx_end-1, idx_start, -1
               if (kdx < idx_end-1) then
                  refval = matrix(kdx+2, kdx+1)
                  zeroval = matrix(kdx+2, kdx)
               end if

               givens = Givens_CS(refval=refval, zeroval=zeroval)
               if (kdx == idx_end-1) then
                  if (all(givens == identity_2d)) then
                     ! If the givens rotation is the same as the identity matrix in the used precision
                     ! then explicitly set the investigated value to zero and continue with the next value
                     matrix(idx_end, idx_end-1) = 0.0_dp
                     exit
                  end if
               end if

               ! Instead of creating the full matrices for a Givens rotation `\mathbf{G}\cdot\mathbf{A}\cdot\mathbf{G}^T`, only use the relevant
               ! slices (only row and column of kdx and kdx+1 change, but two commands seem necessary)
               matrix(:, [kdx, kdx+1]) = matmul(matrix(:, [kdx, kdx+1]), givens)
               matrix([kdx, kdx+1], :) = matmul(transpose(givens), matrix([kdx, kdx+1], :))

               eigenvec_mat(:, [kdx, kdx+1]) = matmul(eigenvec_mat(:, [kdx, kdx+1]), givens)
            end do
         end do
      end do

      ! All nondiagonal elements should be zero or close to it. Explicitly set them to zero and report success
      do idx = 1, nel
         refval = matrix(idx, idx)
         matrix(:, idx) = 0.0_dp
         matrix(idx, idx) = refval
      end do
      success = .True.
   end subroutine QL_Decomposition


   ! --------------------------------------------------------------- !
   subroutine Eigendecomposition(mat, nel, eigenvalues, eigenvectors)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      integer, intent(in) :: nel
      real(dp), dimension(nel, nel), intent(in) :: mat
      real(dp), dimension(nel, nel), intent(out) :: eigenvalues, eigenvectors
      ! ------------------------------------------------------------ !
      real(dp), dimension(nel, nel) :: temp_mat, trans_mat
      logical :: success

      temp_mat = mat
      call Tridiagonal_Matrix(matrix=temp_mat, nel=nel, trans_mat=trans_mat)
      call QL_Decomposition(matrix=temp_mat, nel=nel, eigenvec_mat=eigenvectors, success=success)

      if (.not. success) then
         call Write_Error_And_Exit('Eigendecomposition: Failed QL decomposition of' // char(10) &
            // Formatval('mat', mat))
      end if

      eigenvectors = matmul(trans_mat, eigenvectors)
      eigenvalues = temp_mat
   end subroutine Eigendecomposition


   ! --------------------------------------------------------------- !
   function Matrix_Exponential(mat) result(mat_out)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: mat
      real(dp), dimension(6) :: mat_out
      ! ------------------------------------------------------------ !
      real(dp), dimension(3, 3) :: temp_mat, eigenval_diag, eigenvec
      integer :: idx

      ! `e^{a\mathbf{M}}` is meant to apply the `e^ax`-function to each entry `x` of the diagonal matrix `\mathbf{M}`
      ! If `\mathbf{M}` is not diagonal, decompose its eigenvalues/eigenvectors first, so that `e^{a\mathbf{M}} = \mathbf{V}\cdot\mathrm{diag}\left(e^{a\lambda_1},\dots,e^{a\lambda_n}\right)\mathbf{V}^{-1}`
      if (Is_Diagonal(mat)) then
         mat_out = 0.0_dp
         do idx = 1, 3
            mat_out(idx) = exp(mat(idx))
         end do
      else
         temp_mat = Vec_To_Mat33(vec6=mat)
         call Eigendecomposition(mat=temp_mat, nel=3, eigenvalues=eigenval_diag, eigenvectors=eigenvec)

         temp_mat = 0.0_dp
         do idx = 1, 3
            temp_mat(idx, idx) = exp(eigenval_diag(idx, idx))
         end do
         temp_mat = matmul(eigenvec, matmul(temp_mat, transpose(eigenvec)))
         mat_out = Mat33_To_Vec(mat33=temp_mat)
      end if
   end function Matrix_Exponential


   ! --------------------------------------------------------------- !
   pure function Is_Diagonal(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      logical :: Is_Diagonal
      ! ------------------------------------------------------------ !
      if (any(abs(vec6(4:6)) > setting_epsilon)) then
         Is_Diagonal = .False.
      else
         Is_Diagonal = .True.
      end if
   end function Is_Diagonal


   ! --------------------------------------------------------------- !
   pure function Value_In_Interval(value, limits)
   ! --------------------------------------------------------------- !
      real(dp), intent(in) :: value
      real(dp), dimension(2), intent(in) :: limits
      logical :: Value_In_Interval
      ! ------------------------------------------------------------ !
      Value_In_Interval = .False.
      if ((value >= minval(limits)) .and. (value <= maxval(limits))) then
         Value_In_Interval = .True.
      end if
   end function Value_In_Interval


   ! --------------------------------------------------------------- !
   subroutine Abort_If_Not_In_Interval(name, number, limits)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: number
      real(dp), dimension(2), intent(in) :: limits
      ! ------------------------------------------------------------ !
      if (.not. Value_In_Interval(number, limits)) then
         call Write_Error_And_Exit('Abort_If_Not_In_Interval: ' // &
            Formatval(name, number) // ' not in valid interval')
      end if
   end subroutine Abort_If_Not_In_Interval


   ! --------------------------------------------------------------- !
   pure subroutine Set_Element_In_Tensor(tens, idx, jdx, kdx, ldx, val)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(inout) :: tens
      integer, intent(in) :: idx, jdx, kdx, ldx
      real(dp), intent(in) :: val
      ! ------------------------------------------------------------ !
      tens(ref_elements(idx, jdx), ref_elements(kdx, ldx)) = val
   end subroutine Set_Element_In_Tensor


   ! --------------------------------------------------------------- !
   pure function Get_Matrix_From_Tensorpart(tens, idx, jdx)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(in) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(6) :: Get_Matrix_From_Tensorpart
      ! ------------------------------------------------------------ !
      Get_Matrix_From_Tensorpart = tens(:, ref_elements(idx, jdx))
   end function Get_Matrix_From_Tensorpart


   ! --------------------------------------------------------------- !
   pure subroutine Set_Matrix_In_Tensorpart(tens, idx, jdx, mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6, 6), intent(inout) :: tens
      integer, intent(in) :: idx, jdx
      real(dp), dimension(6), intent(in) :: mat
      ! ------------------------------------------------------------ !
      tens(:, ref_elements(idx, jdx)) = mat
   end subroutine Set_Matrix_In_Tensorpart


   ! --------------------------------------------------------------- !
   pure function Get_Elements_From_Matrixlist(matlist, nel, idx, jdx)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 6), intent(in) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel) :: Get_Elements_From_Matrixlist
      ! ------------------------------------------------------------ !
      Get_Elements_From_Matrixlist = matlist(:, ref_elements(idx, jdx))
   end function Get_Elements_From_Matrixlist


   ! --------------------------------------------------------------- !
   pure subroutine Set_Elements_In_Matrixlist(matlist, nel, idx, jdx, elemlist)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: nel
      real(dp), dimension(nel, 6), intent(inout) :: matlist
      integer, intent(in) :: idx, jdx
      real(dp), dimension(nel), intent(in) :: elemlist
      ! ------------------------------------------------------------ !
      matlist(:, ref_elements(idx, jdx)) = elemlist
   end subroutine Set_Elements_In_Matrixlist


   ! --------------------------------------------------------------- !
   pure subroutine Ref_Index(ref_idx, idx, jdx)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: ref_idx
      integer, intent(out) :: idx, jdx
      ! ------------------------------------------------------------ !
      idx = ref_indices(1, ref_idx)
      jdx = ref_indices(2, ref_idx)
   end subroutine Ref_Index


   ! --------------------------------------------------------------- !
   pure subroutine Perturbate(mat, idx, jdx, theta)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(inout) :: mat
      integer, intent(in) :: idx, jdx
      real(dp), intent(in) :: theta
      ! ------------------------------------------------------------ !
      integer :: comb_idx

      comb_idx = ref_elements(idx, jdx)
      if (idx == jdx) then
         mat(comb_idx) = mat(comb_idx) + theta
      else
         mat(comb_idx) = mat(comb_idx) + 0.5_dp*theta
      end if
   end subroutine Perturbate


   ! --------------------------------------------------------------- !
   pure function Mat33_To_Vec(mat33) result(vec6)
   ! --------------------------------------------------------------- !
      real(dp), dimension(3, 3), intent(in) :: mat33
      real(dp), dimension(6) :: vec6
      ! ------------------------------------------------------------ !
      vec6 = [mat33(1, 1), mat33(2, 2), mat33(3, 3), &
         0.5_dp*(mat33(1, 2) + mat33(2, 1)), 0.5_dp*(mat33(2, 3) + mat33(3, 2)), 0.5_dp*(mat33(1, 3) + mat33(3, 1))]
   end function Mat33_To_Vec


   ! --------------------------------------------------------------- !
   pure function Vec_To_Mat33(vec6) result(mat33)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: vec6
      real(dp), dimension(3, 3) :: mat33
      ! ------------------------------------------------------------ !
      mat33 = reshape([vec6(1), vec6(4), vec6(6), &
                       vec6(4), vec6(2), vec6(5), &
                       vec6(6), vec6(5), vec6(3)], [3, 3])
   end function Vec_To_Mat33


   ! --------------------------------------------------------------- !
   pure function Vec9_To_Mat(vec9) result(mat)
   ! --------------------------------------------------------------- !
      real(dp), dimension(9), intent(in) :: vec9
      real(dp), dimension(6) :: mat
      ! ------------------------------------------------------------ !
      mat = [vec9(1:3), 0.5_dp*(vec9(4)+vec9(5)), 0.5_dp*(vec9(6)+vec9(7)), 0.5_dp*(vec9(8)+vec9(9))]
   end function Vec9_To_Mat


   ! --------------------------------------------------------------- !
   pure function Mat_To_Vec9(mat) result(vec9)
   ! --------------------------------------------------------------- !
      real(dp), dimension(6), intent(in) :: mat
      real(dp), dimension(9) :: vec9
      ! ------------------------------------------------------------ !
      vec9 = [mat(1:4), mat(4), mat(5), mat(5), mat(6), mat(6)]
   end function Mat_To_Vec9


   ! --------------------------------------------------------------- !
   pure function Pack_States(stress, jac_stress, statevariables, jac_statevariables) result(output_states)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, setting_max_internal_states
      !
      real(dp), dimension(6), intent(in) :: stress
      real(dp), dimension(6, 6), intent(in) :: jac_stress
      real(dp), dimension(setting_num_statevariables), intent(in) :: statevariables
      real(dp), dimension(setting_num_statevariables, 6), intent(in) :: jac_statevariables
      real(dp), dimension(setting_max_internal_states) :: output_states
      ! ------------------------------------------------------------ !
      output_states(1:6) = reshape(stress, [6])
      output_states(1+6:(6+1)*6) = reshape(jac_stress, [6*6])
      output_states(1+(6+1)*6:(6+1)*6+setting_num_statevariables) = statevariables
      output_states(1+(6+1)*6+setting_num_statevariables:(6+1)*6+(6+1)*setting_num_statevariables) = &
         reshape(jac_statevariables, [setting_num_statevariables*6])
   end function Pack_States


   ! --------------------------------------------------------------- !
   pure subroutine Unpack_States(input_states, stress, jac_stress, statevariables, jac_statevariables)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, setting_max_internal_states
      !
      real(dp), dimension(setting_max_internal_states), intent(in) :: input_states
      real(dp), dimension(6), intent(out) :: stress
      real(dp), dimension(6, 6), intent(out) :: jac_stress
      real(dp), dimension(setting_num_statevariables), intent(out) :: statevariables
      real(dp), dimension(setting_num_statevariables, 6), intent(out) :: jac_statevariables
      ! ------------------------------------------------------------ !
      stress = reshape(input_states(1:6), [6])
      jac_stress = reshape(input_states(1+6:(6+1)*6), [6, 6])
      statevariables = input_states(1+(6+1)*6:(6+1)*6+setting_num_statevariables)
      jac_statevariables = reshape( &
         input_states(1+(6+1)*6+setting_num_statevariables:(6+1)*6+(6+1)*setting_num_statevariables), &
         [setting_num_statevariables, 6])
   end subroutine Unpack_States


   ! --------------------------------------------------------------- !
   pure function Import_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat has to be adjusted beforehand
      real(dp), dimension(num_dimensions+num_shear), intent(in) :: mat
      real(dp), dimension(6) :: mat_out
      ! ------------------------------------------------------------ !
      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      mat_out(4:(3 + num_shear)) = mat((num_dimensions + 1):(num_dimensions + num_shear))
   end function Import_Matrix


   ! --------------------------------------------------------------- !
   pure function Export_Matrix(mat, num_dimensions, num_shear) result(mat_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of mat_out may need to be adjusted hereafter
      real(dp), dimension(6), intent(in) :: mat
      real(dp), dimension(num_dimensions+num_shear) :: mat_out
      ! ------------------------------------------------------------ !
      mat_out = 0.0_dp
      mat_out(1:num_dimensions) = mat(1:num_dimensions)
      mat_out((num_dimensions + 1):(num_dimensions + num_shear)) = mat(4:(3 + num_shear))
   end function Export_Matrix


   ! --------------------------------------------------------------- !
   pure function Export_Tensor(tens, num_dimensions, num_shear) result(tens_out)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: num_dimensions, num_shear               ! The precision of tens_out may need to be adjusted hereafter
      real(dp), dimension(6, 6), intent(in) :: tens
      real(dp), dimension(num_dimensions+num_shear, num_dimensions+num_shear) :: tens_out
      ! ------------------------------------------------------------ !
      integer :: nel

      nel = num_dimensions + num_shear
      tens_out = 0.0_dp
      tens_out(1:num_dimensions, 1:num_dimensions) = tens(1:num_dimensions, 1:num_dimensions)
      tens_out((num_dimensions + 1):nel, 1:num_dimensions) = tens(4:(3 + num_shear), 1:num_dimensions)
      tens_out(1:num_dimensions, (num_dimensions + 1):nel) = tens(1:num_dimensions, 4:(3 + num_shear))
      tens_out((num_dimensions + 1):nel, (num_dimensions + 1):nel) = tens(4:(3 + num_shear), 4:(3 + num_shear))
   end function Export_Tensor
end module Math_Operations


! ==================================================================================================================== !
module Constitutive_Model_Baseclass
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   implicit none

   private
   public :: Constitutive_Model

   ! --------------------------------------------------------------- !
   type, abstract :: Constitutive_Model
   ! --------------------------------------------------------------- !
      real(dp) :: overall_dt
      real(dp), dimension(6) :: dot_strain
      logical :: calculate_jacobian
      logical :: provide_jacobian
      real(dp), dimension(setting_num_statevariables) :: direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask

      contains

      procedure(initialize_interface), deferred :: Initialize
      procedure(calculate_dot_state_interface), deferred :: Calculate_Dot_State

      procedure :: Base_Initialization
      procedure :: Get_Dot_State
      procedure :: Set_Values
      procedure :: Get_Direct_Variables
      procedure :: Get_Direct_Variables_Mask
      procedure :: Approximate_Jacobian
      procedure :: Consistent_Jacobian
      procedure :: Elasticity
      procedure, nopass :: Is_Valid_Stress_State
      procedure, nopass :: Get_Factor_F
      procedure :: Intergranular_Strain
   end type

   ! --------------------------------------------------------------- !
   abstract interface
   ! --------------------------------------------------------------- !
      subroutine initialize_interface(this, params, calculate_jacobian, firstcall)
         import
         class(Constitutive_Model), intent(inout) :: this
         real(dp), dimension(:), intent(in) :: params
         logical, intent(in) :: calculate_jacobian, firstcall
      end subroutine initialize_interface

      function calculate_dot_state_interface(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
         import
         class(Constitutive_Model), intent(inout) :: this
         real(dp), intent(in) :: ref_dt
         real(dp), intent(in) :: cur_time
         real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
         real(dp), dimension(6), intent(in) :: dot_strain
         real(dp), dimension(setting_max_internal_states) :: dot_state
      end function calculate_dot_state_interface
   end interface


   contains


   ! --------------------------------------------------------------- !
   subroutine Base_Initialization(this, calculate_jacobian, provide_jacobian)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(inout) :: this
      logical, intent(in) :: calculate_jacobian, provide_jacobian
      ! ------------------------------------------------------------ !
      this%calculate_jacobian = calculate_jacobian
      this%provide_jacobian = provide_jacobian
      this%direct_variables = 0.0_dp
      this%direct_variables_mask = .False.
   end subroutine Base_Initialization


   ! --------------------------------------------------------------- !
   function Get_Dot_State(this, ref_dt, cur_time, cur_state) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian, setting_epsilon
      use Math_Operations, only: Norm
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      dot_state = this%Calculate_Dot_State(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_strain=this%dot_strain)

      if ((this%calculate_jacobian) .and. ((setting_numerical_jacobian) .or. (.not. this%provide_jacobian))) then
         ! Calculate the jacobian if requested (and either the numerical calculation is selected or the
         ! constitutive model does not provide the jacobian). Use approximation if strain increment is (almost) zero
         if (Norm(this%dot_strain) < setting_epsilon) then
            call this%Approximate_Jacobian(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_state=dot_state)
         else
            call this%Consistent_Jacobian(ref_dt=ref_dt, cur_time=cur_time, cur_state=cur_state, dot_state=dot_state)
         end if
      end if
   end function Get_Dot_State


   ! --------------------------------------------------------------- !
   subroutine Set_Values(this, overall_dt, dot_strain)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: overall_dt
      real(dp), dimension(6), intent(in) :: dot_strain
      ! ------------------------------------------------------------ !
      this%overall_dt = overall_dt
      this%dot_strain = dot_strain
   end subroutine Set_Values


   ! --------------------------------------------------------------- !
   subroutine Get_Direct_Variables(this, direct_variables)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(in) :: this
      real(dp), dimension(setting_num_statevariables), intent(out) :: direct_variables
      ! ------------------------------------------------------------ !
      direct_variables = this%direct_variables
   end subroutine Get_Direct_Variables


   ! --------------------------------------------------------------- !
   pure function Get_Direct_Variables_Mask(this) result(direct_var_mask)
   ! --------------------------------------------------------------- !
      class(Constitutive_Model), intent(in) :: this
      logical, dimension(setting_num_statevariables) :: direct_var_mask
      ! ------------------------------------------------------------ !
      direct_var_mask = this%direct_variables_mask
   end function Get_Direct_Variables_Mask


   ! --------------------------------------------------------------- !
   subroutine Approximate_Jacobian(this, ref_dt, cur_time, cur_state, dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian_disturbance
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Perturbate, Set_Matrix_In_Tensorpart
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states), intent(inout) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: tmp_stress, dot_stress, dot_stress1_pert, dot_stress2_pert
      real(dp), dimension(6) :: dot_strain_pert, tmp_jac_mat
      real(dp), dimension(setting_num_statevariables) :: tmp_statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: tmp_jacobian, dot_jacobian
      real(dp), dimension(setting_num_statevariables, 6) :: tmp_sta_jacobian, dot_sta_jacobian
      real(dp), dimension(setting_max_internal_states) :: dot_state1_pert, dot_state2_pert
      integer :: idx_loop, idx, jdx

      call Unpack_States(input_states=dot_state, stress=dot_stress, jac_stress=dot_jacobian, &
         statevariables=dot_statevariables, jac_statevariables=dot_sta_jacobian)
      ! The values of dot_jacobian and dot_sta_jacobian are assumed to be initialized as zero in Calculate_Results()

      do idx_loop = 1, 6
         call Ref_Index(ref_idx=idx_loop, idx=idx, jdx=jdx)

         dot_strain_pert = 0.0_dp
         call Perturbate(mat=dot_strain_pert, idx=idx, jdx=jdx, theta=-setting_numerical_jacobian_disturbance)

         dot_state1_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=cur_state, dot_strain=dot_strain_pert)
         call Unpack_States(input_states=dot_state1_pert, stress=dot_stress1_pert, jac_stress=tmp_jacobian, &
            statevariables=tmp_statevariables, jac_statevariables=tmp_sta_jacobian)

         dot_state2_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=cur_state, dot_strain=-dot_strain_pert)
         call Unpack_States(input_states=dot_state2_pert, stress=dot_stress2_pert, jac_stress=tmp_jacobian, &
            statevariables=tmp_statevariables, jac_statevariables=tmp_sta_jacobian)

         tmp_jac_mat= (dot_stress2_pert - dot_stress1_pert) / (2.0_dp*setting_numerical_jacobian_disturbance)
         call Set_Matrix_In_Tensorpart(tens=dot_jacobian, idx=idx, jdx=jdx, mat=tmp_jac_mat)
      end do

      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jacobian, statevariables=dot_statevariables, &
         jac_statevariables=dot_sta_jacobian)
   end subroutine Approximate_Jacobian


   ! --------------------------------------------------------------- !
   subroutine Consistent_Jacobian(this, ref_dt, cur_time, cur_state, dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian_disturbance
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Perturbate, &
                                 Get_Matrix_From_Tensorpart, Get_Elements_From_Matrixlist, &
                                 Set_Matrix_In_Tensorpart, Set_Elements_In_Matrixlist
      !
      class(Constitutive_Model), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(setting_max_internal_states), intent(inout) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, cur_stress_pert, dot_stress_pert
      real(dp), dimension(6) :: dot_strain_pert, tmp_jac_mat
      real(dp), dimension(setting_num_statevariables) :: cur_statevariables, dot_statevariables, tmp_jac_vec, &
                                                         cur_statevariables_pert, dot_statevariables_pert
      real(dp), dimension(6, 6) :: cur_jacobian, dot_jacobian, tmp_jacobian
      real(dp), dimension(setting_num_statevariables, 6) :: cur_sta_jacobian, dot_sta_jacobian, tmp_sta_jacobian
      real(dp), dimension(setting_max_internal_states) :: tmp_state, dot_state_pert
      integer :: idx_loop, idx, jdx

      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=cur_jacobian, &
         statevariables=cur_statevariables, jac_statevariables=cur_sta_jacobian)
      call Unpack_States(input_states=dot_state, stress=dot_stress, jac_stress=dot_jacobian, &
         statevariables=dot_statevariables, jac_statevariables=dot_sta_jacobian)
      ! The values of dot_jacobian and dot_sta_jacobian are assumed to be initialized as zero in the constitutive model(s)

      do idx_loop = 1, 6
         call Ref_Index(ref_idx=idx_loop, idx=idx, jdx=jdx)

         tmp_jac_mat = Get_Matrix_From_Tensorpart(tens=cur_jacobian, idx=idx, jdx=jdx)
         tmp_jac_vec = Get_Elements_From_Matrixlist(matlist=cur_sta_jacobian, nel=setting_num_statevariables, &
            idx=idx, jdx=jdx)

         cur_stress_pert = cur_stress + setting_numerical_jacobian_disturbance * tmp_jac_mat
         cur_statevariables_pert = cur_statevariables + setting_numerical_jacobian_disturbance * tmp_jac_vec
         tmp_state = Pack_States(stress=cur_stress_pert, jac_stress=cur_jacobian, &
            statevariables=cur_statevariables_pert, jac_statevariables=cur_sta_jacobian)

         dot_strain_pert =  this%dot_strain
         call Perturbate(mat=dot_strain_pert, idx=idx, jdx=jdx, theta=setting_numerical_jacobian_disturbance)

         dot_state_pert = this%Calculate_Dot_State(ref_dt=ref_dt, &
            cur_time=cur_time, cur_state=tmp_state, dot_strain=dot_strain_pert)

         call Unpack_States(input_states=dot_state_pert, stress=dot_stress_pert, jac_stress=tmp_jacobian, &
            statevariables=dot_statevariables_pert, jac_statevariables=tmp_sta_jacobian)

         tmp_jac_mat = (dot_stress_pert - dot_stress) / setting_numerical_jacobian_disturbance
         tmp_jac_vec = (dot_statevariables_pert - dot_statevariables) / setting_numerical_jacobian_disturbance

         call Set_Matrix_In_Tensorpart(tens=dot_jacobian, idx=idx, jdx=jdx, mat=tmp_jac_mat)
         call Set_Elements_In_Matrixlist(matlist=dot_sta_jacobian, nel=setting_num_statevariables, &
            idx=idx, jdx=jdx, elemlist=tmp_jac_vec)
      end do

      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jacobian, statevariables=dot_statevariables, &
         jac_statevariables=dot_sta_jacobian)
   end subroutine Consistent_Jacobian


   ! --------------------------------------------------------------- !
   pure subroutine Elasticity(this, youngs_modulus, nu, dot_strain, dot_stress, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_numerical_jacobian
      use Math_Operations, only: Double_Contraction42, const_identity4d_sym, const_identity4d_tr
      !
      class(Constitutive_Model), intent(in) :: this
      real(dp), intent(in) :: youngs_modulus, nu
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(6), intent(out) :: dot_stress
      real(dp), dimension(6, 6), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      associate(Youngsmat => jacobian)
         ! Using the formula is not the most efficient way, but it is independent of the used representation
         Youngsmat = youngs_modulus/(1.0_dp + nu) &                  ! `\mathcal{E} = \frac{E}{1+\nu}\left(\mathcal{I}^\mathrm{sym} + \frac{\nu}{1-2\nu}\mathcal{I}^\mathrm{tr}\right)`
                   * const_identity4d_sym + youngs_modulus*nu/(1.0_dp - nu - 2.0_dp*nu**2)*const_identity4d_tr
         dot_stress = Double_Contraction42(Youngsmat, dot_strain)    ! `\mathbf{\dot{T}} = \mathcal{E}:\mathbf{D}`
      end associate

      if ((.not. this%calculate_jacobian) .or. (setting_numerical_jacobian)) then
         jacobian = 0.0_dp
      end if
   end subroutine Elasticity


   ! --------------------------------------------------------------- !
   pure function Is_Valid_Stress_State(stress) result(in_limits)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_very_small_stress
      use Math_Operations, only: Trace
      !
      real(dp), dimension(6), intent(in) :: stress
      logical :: in_limits
      ! ------------------------------------------------------------ !
      real(dp) :: trT

      ! NOTE: Currently this function checks the trace only. For a more reliable check the individual diagonal components should be checked
      !       and possibly also the actually used limit criterion (Mohr-Coulomb, Matsuoka-Nakai, ...)
      in_limits = .True.
      trT = Trace(stress)
      if (trT > setting_very_small_stress) then
         in_limits = .False.
      end if
   end function Is_Valid_Stress_State


   ! --------------------------------------------------------------- !
   pure function Get_Factor_F(stress_dless_dev)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: const_root2, const_root3, const_root6, Norm, Determinant, Double_Contraction22
      !
      real(dp), dimension(6), intent(in) :: stress_dless_dev
      real(dp) :: Get_Factor_F
      ! ------------------------------------------------------------ !
      real(dp) :: tan_psi, stress_contraction, cos_3theta

      tan_psi = const_root3*Norm(stress_dless_dev)                   ! (2.67) (1) of Niemunis (2003): `\tan{(\psi)} = \sqrt{3}||\hat{\mathbf{T}}^\mathrm{dev}||`
      stress_contraction = Double_Contraction22(stress_dless_dev, stress_dless_dev)

      cos_3theta = 1.0_dp                                            ! Assume term equals one in case denominator `\hat{\mathbf{T}}^\mathrm{dev}:\hat{\mathbf{T}}^\mathrm{dev}` is zero
      if (stress_contraction > setting_epsilon) then
         ! For the given definition of `\hat{\mathbf{T}}^\mathrm{dev}`, the terms `\tr{(\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev})}` and `3\det{(\hat{\mathbf{T}}^\mathrm{dev})}` are equal
         cos_3theta = -const_root6 * 3.0_dp &                        ! (2.67) (2) of Niemunis (2003): `\cos{(3\theta)} = -\sqrt{6}\frac{\tr{(\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev}\cdot\hat{\mathbf{T}}^\mathrm{dev})}}{\left[\hat{\mathbf{T}}^\mathrm{dev}:\hat{\mathbf{T}}^\mathrm{dev}\right]^{1.5}}`
                    * Determinant(stress_dless_dev) / (stress_contraction**1.5_dp)
      end if
      cos_3theta = max(-1.0_dp, min(1.0_dp, cos_3theta))             ! Restrict `\cos{(3\theta)}` in interval `[-1, 1]`

      Get_Factor_F = sqrt(abs(tan_psi*tan_psi/8.0_dp &               ! Modified (2.66) of Niemunis (2003):
                   + (2.0_dp - tan_psi*tan_psi) &                    ! `F = \sqrt{|\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}|} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
                   /(2.0_dp + const_root2*tan_psi*cos_3theta))) - tan_psi/(2.0_dp*const_root2)
   end function Get_Factor_F


   ! --------------------------------------------------------------- !
   pure subroutine Intergranular_Strain(this, L_mat, fdN_mat, D_vis, hypoplastic, igran_strain, &
      R_max, m_T, m_R, beta_R, chi, dt, dot_strain, dot_stress, dot_igran_strain, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, setting_numerical_jacobian
      use Math_Operations, only: Nonzero_Division, Norm, Double_Contraction22, Double_Contraction42, &
                                 Double_Contraction44, Dyadic_Product22, const_identity4d_sym
      !
      class(Constitutive_Model), intent(in) :: this
      real(dp), dimension(6, 6), intent(in) :: L_mat
      real(dp), dimension(6), intent(in) :: fdN_mat, D_vis, igran_strain
      logical, intent(in) :: hypoplastic
      real(dp), intent(in) :: R_max, m_T, m_R, beta_R, chi, dt
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(6), intent(out) :: dot_stress, dot_igran_strain
      real(dp), dimension(6, 6), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      real(dp) :: norm_h, norm_D, rho, load_dir, r_c, ref_dt, sub_dt, sub_dtmin, done_dt
      real(dp), dimension(6) :: h_dir, new_igran_strain
      real(dp), dimension(6, 6) :: L_mod
      logical :: precheck, breakloop

      ref_dt = dt
      ! In case dt is zero, use ref_dt = 1 to prevent division by zero
      if (ref_dt < setting_epsilon) then
         ref_dt = 1.0_dp
      end if

      dot_stress = 0.0_dp
      jacobian = 0.0_dp

      done_dt = 0.0_dp
      new_igran_strain = igran_strain

      norm_D = Norm(dot_strain)
      ! Use given dt (which might be zero) only here to calculate precheck
      precheck = (norm_D*dt <= 0.1_dp*R_max)
      breakloop = .False.

      associate(M_mat => L_mod)
         igranloop: &
         do
            norm_h = Norm(new_igran_strain)
            rho = norm_h/R_max                                       ! (4.5) of Niemunis (2003): `\rho = \frac{||\mathbf{h}||}{R_\mathrm{max}}` with `0 \leq \rho \leq 1`

            load_dir = Double_Contraction22(new_igran_strain, &      ! Proportional to loading direction `\vec{\mathbf{h}}:\mathbf{D}`
                       dot_strain)

            if (precheck) then
               sub_dt = ref_dt
               breakloop = .True.
            else
               ! This branch can't be reached if either dt or norm_D is zero
               sub_dtmin = 0.1_dp*R_max/norm_D

               if (load_dir > 0.0_dp) then
                  if (rho < 0.99_dp) then
                     sub_dt = min(1.01_dp*ref_dt, abs(0.1_dp*R_max/((1.0_dp - rho**beta_R)*norm_D)))
                  else
                     sub_dt = 1.01_dp*ref_dt
                  end if
               else
                  sub_dt = sub_dtmin
               end if

               if ((rho > 0.01_dp) .and. (load_dir >= 0.0_dp)) then
                  ! if norm_h is zero, so is rho
                  sub_dt = sub_dtmin + (sub_dt - sub_dtmin)*load_dir/(norm_h*norm_D)
               end if

               if (done_dt + sub_dt >= ref_dt) then
                  sub_dt = ref_dt - done_dt
                  breakloop = .True.
               end if
            end if

            h_dir = Nonzero_Division(val=new_igran_strain, &         ! Strain direction `\vec{\mathbf{h}} = \frac{\mathbf{h}}{||\mathbf{h}||}`
               fac=norm_h)

            r_c = rho**Chi
            L_mod = Dyadic_Product22(Double_Contraction42(L_mat, &   ! `\mathcal{L}_\mathrm{mod} = \mathcal{L}:\vec{\mathbf{h}}\otimes\vec{\mathbf{h}}`
                    h_dir), h_dir)

            if (load_dir > 0.0_dp) then                              ! if `\vec{\mathbf{h}}:\mathbf{D} > 0` (loading branch) then use
               M_mat = (r_c*m_T + (1.0_dp - r_c)*m_R)*L_mat &        ! upper branch of (4.11) of Niemunis (2003):
                     + r_c*(1.0_dp - m_T)*L_mod                      ! `\mathcal{M} = \left[\rho^\chi m_T + (1 - \rho^\chi)m_R\right]\mathcal{L} + \rho^\chi(1 - m_T)\mathcal{L}_\mathrm{mod}`
               if (hypoplastic) then
                  M_mat = M_mat &                                    ! with hypoplasticity adjust `\mathcal{M} = \mathcal{M} + \rho^\chi \mathbf{N}\otimes \vec{\mathbf{h}}`
                        + r_c*Dyadic_Product22(fdN_mat, h_dir)
               end if

               dot_igran_strain = Double_Contraction42( &            ! upper branch of (4.12) of Niemunis (2003): `\overset{\circ}{\mathbf{h}} = \left(\mathcal{I} - \vec{\mathbf{h}}\otimes\vec{\mathbf{h}} \rho^{\beta_R}\right):\mathbf{D}`
                                  const_identity4d_sym - rho**beta_R*Dyadic_Product22(h_dir, h_dir), dot_strain)
            else
               M_mat = (r_c*m_T + (1.0_dp - r_c)*m_R)*L_mat &        ! lower branch of (4.11) of Niemunis (2003):
                     + r_c*(m_R - m_T)*L_mod                         ! `\mathcal{M} = \left[\rho^\chi m_T + (1 - \rho^\chi)m_R\right]\mathcal{L} + \rho^\chi(m_R - m_T)\mathcal{L}_\mathrm{mod}`
               dot_igran_strain = dot_strain                         ! lower branch of (4.12) of Niemunis (2003): `\overset{\circ}{\mathbf{h}} = \mathbf{D}`
            end if

            ! Estimate the integrated intergranular strain and determine its norm. Adjust it, if it is greater than R_max
            new_igran_strain = new_igran_strain + dot_igran_strain*sub_dt
            norm_h = Norm(new_igran_strain)

            if (norm_h > R_max) then
               new_igran_strain = new_igran_strain * R_max/norm_h
            end if

            dot_stress = dot_stress &                                ! Incremental variant of (4.6) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{M}:\mathbf{D}`
                       + Double_Contraction42(M_mat, dot_strain)*sub_dt/ref_dt
            if (.not. hypoplastic) then
               dot_stress = dot_stress &                             ! Incremental variant of (4.121) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{M}:\mathbf{D} - \mathcal{L}:\mathbf{D}^\mathrm{vis}`
                          - Double_Contraction42(L_mat, D_vis)*sub_dt/ref_dt
            end if

            done_dt = done_dt + sub_dt

            if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
               jacobian = jacobian + M_mat*sub_dt/ref_dt
            end if

            if (breakloop) then
               exit igranloop
            end if
         end do igranloop
         dot_igran_strain = (new_igran_strain - igran_strain)/ref_dt
      end associate
   end subroutine Intergranular_Strain
end module Constitutive_Model_Baseclass


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Elasticity_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Elasticity
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_youngs_modulus, param_nu

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval
      !
      class(Elasticity), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)

      ! --- Parameters and description
      this%param_youngs_modulus = params(1)                          ! Young's modulus `E` (in kPa)
      this%param_nu             = params(2)                          ! Poisson's ratio `\nu`

      if (firstcall) then                                            ! Upper limit on Young's modulus is chosen arbitrarily
         call Abort_If_Not_In_Interval('E', this%param_youngs_modulus, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
      end if
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: Pack_States, Unpack_States
      !
      class(Elasticity), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables

      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)

      call this%Elasticity(youngs_modulus=this%param_youngs_modulus, nu=this%param_nu, &
         dot_strain=dot_strain, dot_stress=dot_stress, jacobian=dot_jac_stress)

      ! --- Packing all dot_states in a long vector
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Elasticity_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Hypoplasticity_Wu92_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Hypoplasticity_Wu92
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_C1, param_C2, param_C3, param_C4

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Hypoplasticity_Wu92), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)

      ! --- Parameters
      this%param_C1 = params(1)
      this%param_C2 = params(2)
      this%param_C3 = params(3)
      this%param_C4 = params(4)
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Nonzero_Division, Norm, Trace, Squared, Deviatoric_Part, Dyadic_Product22, &
                                 Double_Contraction42, const_identity4d_sym, Pack_States, Unpack_States
      !
      class(Hypoplasticity_Wu92), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, N_mat, D_dir
      real(dp), dimension(6, 6) :: L_mat
      real(dp) :: trT, cur_voidratio, dot_voidratio
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! change of void ratio `\dot{e} = \left(1+e_i\right)\tr{(\mathbf{D})}`

      trT = Trace(cur_stress)
      if (this%Is_Valid_Stress_State(stress=cur_stress) .and. (abs(trT) > setting_epsilon)) then
         L_mat = this%param_C1*trT*const_identity4d_sym &            ! `\mathcal{L} = C_1 \tr{(\mathbf{T})}\mathcal{I}^{id} + \frac{C_2}{\tr{(\mathbf{T})}}\mathbf{T}\otimes\mathbf{T}`
               + this%param_C2/trT*Dyadic_Product22(cur_stress, cur_stress)
         N_mat = (this%param_C3*Squared(cur_stress) &                ! `\mathbf{N} = \frac{C_3\mathbf{T}^2 + C_4\mathbf{\hat{T}}^2}{\tr{(\mathbf{T})}}`
               + this%param_C4*Squared(Deviatoric_Part(cur_stress)))/trT
         dot_stress = Double_Contraction42(L_mat, dot_strain) &      ! `\overset{\circ}{\mathbf{T}} = \mathcal{L}:\mathbf{D} + \mathbf{N}||\mathbf{D}||`
                    + N_mat*Norm(dot_strain)

         if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
            D_dir = Nonzero_Division(val=dot_strain, fac=Norm(dot_strain))
            dot_jac_stress = L_mat + Dyadic_Product22(N_mat, D_dir)
         end if
      else
         dot_stress = 0.0_dp
      end if

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Hypoplasticity_Wu92_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Hypoplasticity_VW96_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Hypoplasticity_VW96
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_nu, param_h_s, param_n_H, param_e_d0, param_e_c0, param_e_i0, &
                  param_alpha_H, param_beta_H, param_m_T, param_m_R, param_R_max, param_beta_R, param_chi, &
                  param_e_ini
      real(dp) :: derived_a, derived_a2, derived_f_b_lastterm, derived_b2

      contains

      procedure :: Initialize
      procedure :: Valid_State
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Hypoplasticity_VW96), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)
      this%direct_variables_mask(11) = .True.

      ! --- Standard hypoplastic parameters
      this%param_phi_c   = params(1)                                 ! Friction angle `\varphi_c` (in radians)
      this%param_nu      = params(2)                                 ! Poisson ratio `\nu` for increased shear stiffness
      this%param_h_s     = params(3)                                 ! Granulate hardness `h_s` (in kPa)
      this%param_n_H     = params(4)                                 ! Exponent `n_H`
      this%param_e_d0    = params(5)                                 ! Densest void ratio for zero stress `e_{d0}`
      this%param_e_c0    = params(6)                                 ! Critical void ratio for zero stress `e_{c0}`
      this%param_e_i0    = params(7)                                 ! Loosest void ratio for zero stress `e_{i0}`
      this%param_alpha_H = params(8)                                 ! Pycnotropic exponent `\alpha_H`
      this%param_beta_H  = params(9)                                 ! Barotropic exponent `\beta_H`

      ! --- Extended parameters for intergranular strain
      this%param_m_T    = params(10)                                 ! Orthogonal loading history multiplier `m_T`
      this%param_m_R    = params(11)                                 ! Reverse loading history multiplier `m_R`
      this%param_R_max  = params(12)                                 ! Elastic strain range `R_\mathrm{max}`
      this%param_beta_R = params(13)                                 ! Parameter to control evolution of intergranular strain `\beta_R`
      this%param_Chi    = params(14)                                 ! Parameter to control degradation of stiffness `\chi`

      ! --- Additional parameters
      ! param(15) empty due to params structure compatibility with some other routines
      this%param_e_ini  = params(16)                                 ! Void ratio at K0-state `e_\mathrm{ini}`


      ! --- Check valid parameter range
      if (firstcall) then
         ! Limits taken/adjusted from user routine of Niemunis
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
         call Abort_If_Not_In_Interval('h_s', this%param_h_s, [100.0_dp, 1.0e9_dp])
         call Abort_If_Not_In_Interval('n_H', this%param_n_H, [0.0_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('e_d0', this%param_e_d0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [this%param_e_d0, 20.0_dp])
         call Abort_If_Not_In_Interval('e_i0', this%param_e_i0, [this%param_e_c0, 20.0_dp])
         call Abort_If_Not_In_Interval('m_T', this%param_m_T, [0.4_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('m_R', this%param_m_R, [0.4_dp, 20.0_dp])
         call Abort_If_Not_In_Interval('R_max', this%param_R_max, [setting_epsilon, 0.1_dp])
         call Abort_If_Not_In_Interval('beta_R', this%param_beta_R, [setting_epsilon, 100.0_dp])
         call Abort_If_Not_In_Interval('Chi', this%param_Chi, [0.0_dp, 100.0_dp])
      end if

      ! --- Derived parameters
      this%derived_a = const_root3*(3.0_dp - sin(this%param_phi_c)) &! Precalculate (2.65) of Niemunis (2003): `a = \frac{\sqrt{3}\left(3-\sin{(\varphi_c)}\right)}{2\sqrt{2}\sin{(\varphi_c)}}`
                     / (2.0_dp*const_root2*sin(this%param_phi_c))
      this%derived_a2 = this%derived_a**2                            ! `a^2`

      this%derived_f_b_lastterm = 3.0_dp + this%derived_a2 &         ! Precalculate `f_{b,\mathrm{lastterm}} = 3+a^2-a\sqrt{3}\left(\frac{e_{i0}-e_{d0}}{e_{c0}-e_{d0}}\right)^{\alpha}`
                                - this%derived_a*const_root3 &       ! `f_{b,\mathrm{lastterm}}` is the last term of (2.72) in Niemunis (2003)
                                * (((this%param_e_i0 - this%param_e_d0) &
                                / (this%param_e_c0 - this%param_e_d0))**this%param_alpha_H)
      if (this%derived_f_b_lastterm <= 0.0_dp) then
         call Write_Error_And_Exit('Param_Hypoplasticity: Last term of f_b is not greater than zero ' // &
            '(sign of stresses would change)' // char(10) // Formatval('f_b_lastterm', this%derived_f_b_lastterm))
      end if

      this%derived_b2 = (1.0_dp + this%derived_a2/3.0_dp &           ! From (4.183) of Niemunis (2003): `b^2 = \left(1+\frac{a^2}{3} + \frac{a}{\sqrt{3}}\right)\frac{1-2\nu}{1+\nu} - 1`
                      + this%derived_a/const_root3) &                ! (needed for increased shear stiffness)
                      * (1.0_dp - 2.0_dp*this%param_nu)/(1.0_dp + this%param_nu) - 1.0_dp
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Valid_State(this, cur_stress, cur_voidratio, e_d, e_c, e_i, dt, dot_voidratio)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan, Write_Warning, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Nonzero_Division, Trace, Value_In_Interval
      !
      class(Hypoplasticity_VW96), intent(in) :: this
      real(dp), dimension(6), intent(in) :: cur_stress
      real(dp), intent(inout) :: cur_voidratio
      real(dp), intent(out) :: e_d, e_c, e_i
      real(dp), intent(in) :: dt
      real(dp), intent(out) :: dot_voidratio
      logical :: Valid_State
      ! ------------------------------------------------------------ !
      real(dp) :: trT, Bauer_Term, temp_voidratio

      Valid_State = .False.
      dot_voidratio = 0.0_dp
      temp_voidratio = cur_voidratio

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      if (Is_Nan(trT)) then
         call Write_Error_And_Exit('Valid_State: trT is NaN')
      end if

      if (this%Is_Valid_Stress_State(cur_stress)) then
         ! Compute the current void ratios according to Bauer (1996)
         Bauer_Term = exp(-(-trT/this%param_h_s)**this%param_n_H)    ! Exponential compression law by Bauer (1996): `B(\tr{(\mathbf{T})}) = \exp{}[-(\frac{\tr{(\mathbf{T})}}{h_s})^n_H]`
         e_d = this%param_e_d0 * Bauer_Term                          ! Current densest void ratio  `e_d = e_{d0} B(\tr{(\mathbf{T})})`
         e_c = this%param_e_c0 * Bauer_Term                          ! Current critical void ratio `e_c = e_{c0} B(\tr{(\mathbf{T})})`
         e_i = this%param_e_i0 * Bauer_Term                          ! Current loosest void ratio  `e_i = e_{i0} B(\tr{(\mathbf{T})})`

         ! Assuming void ratio is not specified as statevariable try materialparameter
         if (cur_voidratio < setting_epsilon) then
            cur_voidratio = this%param_e_ini * Bauer_Term
         end if

         if (.not. Value_In_Interval(cur_voidratio, [e_d, e_i])) then
            ! Only give feedback if it is considerably off
            if (.not. Value_In_Interval(cur_voidratio, [0.96_dp*e_d, 1.04_dp*e_i])) then
               call Write_Warning('Valid_State: ' // Formatval('void ratio:', cur_voidratio) // &
                  ' out of interval [' // Formatval('e_d:', e_d) // ', ' // Formatval('e_i:', e_i) // ']')
            end if

            cur_voidratio = max(e_d, min(e_i, cur_voidratio))
         end if
         Valid_State = .True.
      end if

      dot_voidratio = Nonzero_Division(val=(cur_voidratio - temp_voidratio), fac=dt)
   end function Valid_State


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, setting_hypo_increased_stiffness, setting_hypo_consistent_f_d, &
                                  Write_Error_And_Exit, setting_min_youngs_modulus
      use Math_Operations, only: const_root3, const_identity4d_sym, const_identity4d_tr, &
                                 Nonzero_Division, Norm, Trace, Dimensionless, Deviatoric_Part, Inverse_Tensor, &
                                 Tensor_Partialtrace, Dyadic_Product22, Double_Contraction22, Double_Contraction42, &
                                 Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Hypoplasticity_VW96), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, T_dless, T_dless_dev
      real(dp) :: cur_voidratio, dot_voidratio, trT, e_d, e_c, e_i, norm_D, cur_param_young
      real(dp) :: m_dir_voidratio, m_dir_stress, factor_F, factor_F2, f_LN, f_e, f_b, f_d_bar, f_d_sign, f_d, yout
      real(dp), dimension(6) :: fdN_mat, igran_strain, dot_igran_strain, D_vis
      real(dp), dimension(6, 6) :: L_mat, L_inv, jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      logical :: success

      dot_igran_strain = 0.0_dp

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      igran_strain = Vec9_To_Mat(vec9=statevariables(2:10))
      cur_param_young = statevariables(11)

      ! NOTE: Valid_State has to return .False. if either cur_voidratio or cur_stress are zero
      !       because they are used in the denominator of some hypoplastic expressions below.
      !       Also dot_voidratio will be initialised in Valid_State. It will be zero if the
      !       current void_ratio is already in valid bounds or non-zero otherwise.
      if (.not. this%Valid_State(cur_stress=cur_stress, cur_voidratio=cur_voidratio, &
         e_d=e_d, e_c=e_c, e_i=e_i, dt=ref_dt, dot_voidratio=dot_voidratio)) then

         if (cur_param_young < setting_min_youngs_modulus) then
            cur_param_young = setting_min_youngs_modulus
         end if

         call this%Elasticity(youngs_modulus=cur_param_young, nu=this%param_nu, dot_strain=dot_strain, &
            dot_stress=dot_stress, jacobian=dot_jac_stress)

         dot_igran_strain = Nonzero_Division(val=-igran_strain, &    ! Reset intergranular strain for replacement model
            fac=ref_dt)

         this%direct_variables(11) = cur_param_young

         ! --- Packing all dot_states in a long vector
         dot_statevariables(1) = dot_voidratio
         dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
         dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
            jac_statevariables=dot_jac_statevariables)
         return
      end if

      dot_voidratio = dot_voidratio + &                              ! change of void ratio `\dot{e} = \left(1+e\right)\tr{(\mathbf{D})}`
                      (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! additionally take into account if `e` was not in valid bounds and
                                                                     ! add the necessary change `\dot{e}_{corr}`

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      T_dless = Dimensionless(cur_stress)                            ! Dimensionless stress tensor `\hat{\mathbf{T}}`
      T_dless_dev = Deviatoric_Part(T_dless)                         ! Deviatoric part of dimensionless stress tensor `\hat{\mathbf{T}}^\mathrm{dev}`

      norm_D = Norm(dot_strain)                                      ! Norm of strain increments `||\mathbf{D}||`

      ! Compute scalar factors `f_e`, `f_b` (both barotropy) and `f_d` (pycnotropy)
      f_e = (e_c/cur_voidratio)**this%param_beta_H                   ! (2.70) of Niemunis (2003): `f_e(\tr{(\mathbf{T})},e) = \left(\frac{e_c}{e}\right)^\beta_H`
      f_b = (this%param_e_i0/this%param_e_c0)**this%param_beta_H &   ! (2.72) of Niemunis (2003): `f_b(\tr{(\mathbf{T})}) = \left(\frac{e_{i0}}{e_{c0}}\right)^\beta_H \frac{h_s}{n_H} \frac{1+e_i}{e_i} \left(\frac{-\tr{(\mathbf{T})}}{h_s}\right)^{1-n_H} \left[f_{b,\mathrm{lastterm}}\right]^{-1}`
          * this%param_h_s/this%param_n_H * (1.0_dp + e_i)/e_i &     ! Therefore `f_b(\tr{(\mathbf{T})}) = \left(\frac{e_{i0}}{e_{c0}}\right)^\beta_H \frac{h_s}{n_H} \frac{1+e_i}{e_i} \left(\frac{-\tr{(\mathbf{T})}}{h_s}\right)^{1-n_H} \left[3+a^2-a\sqrt{3}\left(\frac{e_{i0}-e_{d0}}{e_{c0}-e_{d0}}\right)^{\alpha_H}\right]^{-1}`
          * (-trT/this%param_h_s)**(1.0_dp - this%param_n_H) / this%derived_f_b_lastterm

      ! Either use standard definition of `f_d` as in (2.71) or enforce a consistent lower bound as in (4.223) of Niemunis (2003)
      consistent_f_d : &
      if ((setting_hypo_consistent_f_d) &                            ! The changes of the lower bound would only be applied for `e < e_c`
         .and. (cur_voidratio < e_c)) then

         m_dir_voidratio = 1.0_dp                                    ! Direction of void ratio `M_e^{(d)}` and
         m_dir_stress = -e_d/this%param_h_s * this%param_n_H &       ! direction of stresses `M_T^{(d)}` as in (4.211) of Niemunis (2003)
                      * (-trT/this%param_h_s)**(this%param_n_H-1.0_dp)
         f_d_bar = -(m_dir_voidratio*const_root3 &                   ! (4.221) of Niemunis (2003): `\bar{f}_d = -\frac{M_e^{(d)}\sqrt{3}(1+e)+M_T^{(d)}f_b f_e \frac{3}{\sqrt{3}}(3+a^2)}{M_T^{(d)}f_b f_e 3a}`
                 * (1.0_dp + cur_voidratio) + m_dir_stress*f_b*f_e*const_root3 &
                 * (3.0_dp + this%derived_a2)) / (m_dir_stress*f_b*f_e*3.0_dp*this%derived_a)
         f_d_sign = sign(1.0_dp, cur_voidratio - e_d) &              ! Factor from (4.222) of Niemunis (2003)
                  * (abs(cur_voidratio - e_d) &                      ! `f_{d,\mathrm{sign}} = \left(\frac{e-e_d}{e_c-e_d}\right)^{\alpha_H}` if `e > e_d` and `-\left(\frac{|e-e_d|}{e_c-e_d}\right)^{\alpha_H}` otherwise
                  / (e_c - e_d))**this%param_alpha_H
         f_d = f_d_sign + (1.0_dp - f_d_sign)**5 * f_d_bar           ! (4.223) of Niemunis (2003): `f_d = f_{d,\mathrm{sign}} + \left[1 - f_{d,\mathrm{sign}}\right]^z\bar{f}_d` with `z=5`
      else
         f_d = ((cur_voidratio - e_d) / (e_c - e_d)) &               ! (2.71) of Niemunis (2003): `f_d(\tr{(\mathbf{T})},e) = \left(\frac{e-e_d}{e_c-e_d}\right)^{\alpha_H} = r_e^{\alpha_H}`
             **this%param_alpha_H
      end if consistent_f_d

      ! Compute tensor `\mathcal{L}` and matrix `\mathbf{N}`
      factor_F = this%Get_Factor_F(stress_dless_dev=T_dless_dev)     ! (2.66) of Niemunis (2003): `F=\sqrt{\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
      factor_F2 = factor_F**2
      f_LN = f_b*f_e/Double_Contraction22(T_dless, T_dless)          ! `f_\mathrm{LN} = \frac{f_b f_e}{\hat{\mathbf{T}}:\hat{\mathbf{T}}}`

      L_mat = f_LN*(factor_F2*const_identity4d_sym &                 ! (2.63) of Niemunis (2003): `\mathcal{L} = f_\mathrm{LN}a^2 \cdot \left(\left(\frac{F}{a}\right)^2 \mathcal{I} + \hat{\mathbf{T}}\otimes\hat{\mathbf{T}}\right)`
            + this%derived_a2*Dyadic_Product22(T_dless, T_dless))
      fdN_mat = f_d*f_LN*factor_F*this%derived_a &                   ! Adjusted (2.64) of Niemunis (2003): `f_d\mathbf{N} = f_d f_\mathrm{LN}a^2 \cdot \left(\frac{F}{a}\right)\left(\hat{\mathbf{T}} + \hat{\mathbf{T}}^\mathrm{dev}\right)`
              * (T_dless + T_dless_dev)

      if (setting_hypo_increased_stiffness) then                     ! Increase shear stiffness
         L_mat = L_mat + f_LN*this%derived_b2*const_identity4d_sym & ! Modification of (2.63) of Niemunis (2003):
               - this%derived_b2/3.0_dp*const_identity4d_tr          ! `\mathcal{L}_\mathrm{incr} = f_\mathrm{LN} \cdot \left(\left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}\right)`
         L_inv = (const_identity4d_sym - Dyadic_Product22(T_dless, & ! Adaption of (3.33) of Niemunis (2003):
                 T_dless)/(factor_F2/this%derived_a2 &               ! `\mathcal{L}^{-1} = \frac{1}{f_\mathrm{LN}F^2}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right)`
               + Double_Contraction22(T_dless, T_dless)))/(f_LN*factor_F2)
         fdN_mat = Double_Contraction42(L_mat, &                     ! Forge `f_d\mathbf{N}_\mathrm{incr} = f_d\mathcal{L}_\mathrm{incr}:\left(\mathcal{L}^{-1}:\mathbf{N}\right)`
                   Double_Contraction42(L_inv, fdN_mat))             ! to preserve flow direction
      end if

      if (this%param_m_R <= 2.0_dp) then
         ! ------ Hypoplasticity without intergranular strain ------ !
         dot_stress = Double_Contraction42(L_mat, dot_strain) &      ! (2.61) of Niemunis (2003): `\overset{\circ}{\mathbf{T}} = \mathcal{L}:\mathbf{D} + f_d \mathbf{N}||\mathbf{D}||`
                    + fdN_mat*norm_D
         if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
            if ((setting_hypo_increased_stiffness) .and. (cur_time > setting_epsilon) &
               .and. (norm_D > setting_epsilon)) then

               call Inverse_Tensor(tensor=L_inv, inv_tensor=L_mat, & ! Inverse of `\mathcal{L}_\mathrm{incr} = f_\mathrm{LN} \cdot \left(\left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}\right)`
                  success=success)
               if (.not. success) then
                  call Write_Error_And_Exit('Hypoplasticity: Inversion of L_incr failed')
               end if

               yout = min(0.95_dp, 0.95_dp/(Norm(Double_Contraction42(L_inv, fdN_mat)) + setting_epsilon))
               dot_jac_stress = L_mat + yout*Dyadic_Product22(fdN_mat, dot_strain/norm_D)
            else
               dot_jac_stress = L_mat
            end if
         end if
      else
         ! -------- Hypoplasticity with intergranular strain ------- !
         D_vis = 0.0_dp
         call this%Intergranular_Strain(L_mat=L_mat, fdN_mat=fdN_mat, D_vis=D_vis, hypoplastic=.True., &
            igran_strain=igran_strain, R_max=this%param_R_max, m_T=this%param_m_T, m_R=this%param_m_R, &
            beta_R=this%param_beta_R, chi=this%param_chi, dt=ref_dt, dot_strain=dot_strain, &
            dot_stress=dot_stress, dot_igran_strain=dot_igran_strain, jacobian=dot_jac_stress)
      end if

      ! Estimation of current stiffness for replacement model
      if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
         cur_param_young = Tensor_Partialtrace(dot_jac_stress)/3.0_dp
      else
         cur_param_young = Tensor_Partialtrace(L_mat)/3.0_dp
      end if

      this%direct_variables(11) = cur_param_young

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Hypoplasticity_VW96_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Viscohypoplasticity_Ni03_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states, setting_numerical_jacobian
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Viscohypoplasticity_Ni03
   ! --------------------------------------------------------------- !
      private
      ! NOTE: Currently it is confusing to have four very similar named but different parameters:
      !         - param_beta_b (viscohypoplastic parameter)
      !         - param_beta_R (extended hypoplastic parameter)
      !         - derived_beta_b (internal parameter)
      !         - derived_b2 (internal parameter)

      real(dp) :: param_e_100, param_nu, param_lambda, param_kappa, param_beta_b, param_I_v, param_D_r, &
                  param_phi_c, param_m_T, param_m_R, param_R_max, param_beta_R, param_Chi, param_OCR
      real(dp) :: derived_a, derived_a2, derived_pq_incl, derived_beta_b, derived_b2

      contains

      procedure :: Initialize
      procedure :: Valid_State
      procedure :: Get_Pressure_Preconsol
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Error_And_Exit
      use Debug, only: Formatval
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Viscohypoplasticity_Ni03), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      real(dp) :: K_0

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.True.)
      this%direct_variables_mask(11) = .True.

      ! --- Standard viscohypoplastic parameters
      this%param_e_100  = params(1)                                  ! Void ratio `e_{100}` for `p_{e_{100}}=100 kPa` and chosen `D_r`
      this%param_nu     = params(2)                                  ! Poisson ratio `\nu` for increased shear stiffness
      this%param_lambda = params(3)                                  ! Slope of normal consolidated line in void ratio-pressure diagram `\lambda`
      this%param_kappa  = params(4)                                  ! Slope of reloading branch in void ratio-pressure diagram `\kappa`
      this%param_beta_b = params(5)                                  ! Factor `\beta_b` to describe loading shape
      this%param_I_v    = params(6)                                  ! Viscosity index `I_v`
      this%param_D_r    = params(7)                                  ! Reference creep rate `D_r`
      this%param_phi_c  = params(8)                                  ! Friction angle `\varphi_c` (in radians)

      ! --- Extended parameters for intergranular strain
      this%param_m_T    = params(10)                                 ! Orthogonal loading history multiplier `m_T`
      this%param_m_R    = params(11)                                 ! Reverse loading history multiplier `m_R`
      this%param_R_max  = params(12)                                 ! Elastic strain range `R_\mathrm{max}`
      this%param_beta_R = params(13)                                 ! Parameter to control evolution of intergranular strain `\beta_R`
      this%param_chi    = params(14)                                 ! Parameter to control degradation of stiffness `\chi`

      ! --- Additional viscohypoplastic parameters
      this%param_OCR    = params(15)                                 ! Overconsolidation ratio `OCR` (currently not used)

      ! --- Check valid parameter range
      if (firstcall) then
         ! Limits taken/adjusted from user routine of Niemunis
         call Abort_If_Not_In_Interval('e_100', this%param_e_100, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('nu', this%param_nu, [0.0_dp, 0.5_dp])
         call Abort_If_Not_In_Interval('lambda', this%param_lambda, [0.0_dp, 2.0_dp])
         call Abort_If_Not_In_Interval('kappa', this%param_kappa, [0.0_dp, this%param_lambda])
         call Abort_If_Not_In_Interval('beta_b', this%param_beta_b, [0.0_dp, 1.0_dp])
         call Abort_If_Not_In_Interval('I_v', this%param_I_v, [0.0_dp, 1.0_dp])
         call Abort_If_Not_In_Interval('D_r', this%param_D_r, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('m_T', this%param_m_T, [0.4_dp, 10.0_dp])
         call Abort_If_Not_In_Interval('m_R', this%param_m_R, [0.4_dp, 20.0_dp])
         call Abort_If_Not_In_Interval('R_max', this%param_R_max, [setting_epsilon, 0.1_dp])
         call Abort_If_Not_In_Interval('beta_R', this%param_beta_R, [setting_epsilon, 100.0_dp])
         call Abort_If_Not_In_Interval('chi', this%param_chi, [0.0_dp, 100.0_dp])
      end if

      ! --- Derived parameters
      this%derived_a = const_root3*(3.0_dp - sin(this%param_phi_c)) &! Precalculate (2.65) of Niemunis (2003): `a = \frac{\sqrt{3}\left(3-\sin{(\varphi_c)}\right)}{2\sqrt{2}\sin{(\varphi_c)}}`
                     / (2.0_dp*const_root2*sin(this%param_phi_c))
      this%derived_a2 = this%derived_a**2                            ! `a^2`

      this%derived_pq_incl = 6.0_dp*sin(this%param_phi_c) &          ! inclination in p-q diagram `\frac{6\sin{(\varphi_c)}}{3-\sin{(\varphi_c)}}`
                           / (3.0_dp - sin(this%param_phi_c))

      K_0 = (-2.0_dp - this%derived_a2 + sqrt(36.0_dp &              ! (4.94) of Niemunis (2003): `K_0^{up} = K_0 = \frac{-2-a^2+\sqrt{36+36a^2+a^4}}{16}`
          + 36.0_dp*this%derived_a2 + this%derived_a2**2))/16.0_dp   ! which resembles the upper limit of stress ratio `K_0`
      this%derived_beta_b = 1.0_dp/(this%param_kappa * (1.0_dp &     ! (4.72) of Niemunis (2003): `\beta_b = \frac{1}{\kappa_0\left(1 + a^2/(1 + 2K_0)\right)}`
                          + this%derived_a2 / (1.0_dp + 2.0_dp*K_0)))

      this%derived_b2 = (1.0_dp + this%derived_a2/3.0_dp) &          ! Modified (4.183) of Niemunis (2003): `b^2 = \left(1+\frac{a^2}{3}\right)\frac{1-2\nu}{1+\nu} - 1`
                      * (1.0_dp - 2.0_dp*this%param_nu) &            ! (needed for increased shear stiffness)
                      / (1.0_dp + this%param_nu) - 1.0_dp
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Valid_State(this, cur_stress, cur_voidratio, pressure_equiv)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Is_Nan, Write_Error_And_Exit
      use Math_Operations, only: Trace
      !
      class(Viscohypoplasticity_Ni03), intent(in) :: this
      real(dp), dimension(6), intent(in) :: cur_stress
      real(dp), intent(inout) :: cur_voidratio
      real(dp), intent(out) :: pressure_equiv
      logical :: Valid_State
      ! ------------------------------------------------------------ !
      real(dp) :: trT

      Valid_State = .False.

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      if (Is_Nan(trT)) then
         call Write_Error_And_Exit('Valid_State: trT is NaN')
      end if

      if (cur_voidratio < setting_epsilon) then
         call Write_Error_And_Exit('Valid_State: void ratio invalid (too small) and no default value available')
      end if

      if (this%Is_Valid_Stress_State(cur_stress)) then
         pressure_equiv = 100.0_dp*((1.0_dp + this%param_e_100) &    ! Calculate transformed (4.44) of Niemunis (2003): `T_e = T_{e0}\left(\frac{1+e}{1+e_{100}}\right)^{-\frac{1}{\lambda}}`
                        / (1.0_dp + cur_voidratio)) &                ! with `T_{e0} = 100`
                        **(1.0_dp/this%param_lambda)
         Valid_State = .True.
      end if
   end function Valid_State


   ! --------------------------------------------------------------- !
   function Get_Pressure_Preconsol(this, overcrit, trT)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Debug, only: Formatval
      !
      class(Viscohypoplasticity_Ni03), intent(in) :: this
      real(dp), intent(in) :: overcrit, trT
      real(dp) :: Get_Pressure_Preconsol
      ! ------------------------------------------------------------ !
      real(dp) :: p_mean

      p_mean = -trT/3.0_dp
      if (overcrit > 1.0_dp) then
         Get_Pressure_Preconsol = p_mean/2.0_dp*(1.0_dp + overcrit)*(1.0_dp + this%param_beta_b)
      else
         Get_Pressure_Preconsol = p_mean*(1.0_dp - this%param_beta_b &
                                * sqrt(1.0_dp + overcrit*(this%param_beta_b**2 - 1.0_dp)))/(1.0_dp - this%param_beta_b)
      end if

      if (Get_Pressure_Preconsol < 0.0_dp) then
         call Write_Error_And_Exit('Get_Pressure_Preconsol: Negative preconsolidation pressure ' // &
            Formatval('p_e^+: ', Get_Pressure_Preconsol))
      end if
   end function Get_Pressure_Preconsol


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_min_youngs_modulus, Write_Error_And_Exit
      use Math_Operations, only: const_identity2d, const_identity4d_sym, const_identity4d_tr, &
                                 Nonzero_Division, Trace, Dimensionless, Deviatoric_Part, Norm, Inverse_Tensor, &
                                 Tensor_Partialtrace, Double_Contraction42, Double_Contraction22, Dyadic_Product22, &
                                 Double_Contraction44, Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Viscohypoplasticity_Ni03), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, T_dev, T_dless, T_dless_dev, B_mat, B_dir, D_vis, &
                                dpressure_preconsol_ds, N_mat, igran_strain, dot_igran_strain, N_hat
      real(dp), dimension(6, 6) :: L_hat, L_mat, L_hat_inv, A_term, B_term, C_mat, C_inv, &
                                   K_firstinv, K_first, jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: trT, cur_voidratio, dot_voidratio, factor_F, factor_F2, state_OCR, &
                  cur_param_young, p_mean, q_mean, crit_stress_ratio, eta, eta2, &
                  dtth, pressure_equiv, pressure_preconsol, dpressure_preconsol_dp, dpressure_preconsol_dq
      logical :: success

      ! ATTENTION: This routine is not fully tested and the current implementation will probably not work as intended.
      !            Do not use it unless you have checked it thouroughly and improved it first!

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      igran_strain = Vec9_To_Mat(vec9=statevariables(2:10))
      cur_param_young = statevariables(11)

      if (.not. this%Valid_State(cur_stress=cur_stress, cur_voidratio=cur_voidratio, &
         pressure_equiv=pressure_equiv)) then

         if (cur_param_young < setting_min_youngs_modulus) then
            cur_param_young = setting_min_youngs_modulus
         end if

         call this%Elasticity(youngs_modulus=cur_param_young, nu=this%param_nu, dot_strain=dot_strain, &
            dot_stress=dot_stress, jacobian=dot_jac_stress)

         dot_voidratio = 0.0_dp                                      ! Keep current void ratio
         dot_igran_strain = Nonzero_Division(val=-igran_strain, &    ! Reset intergranular strain for replacement model
            fac=ref_dt)

         this%direct_variables(11) = cur_param_young

         ! --- Packing all dot_states in a long vector
         dot_statevariables(1) = dot_voidratio
         dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
         dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
            jac_statevariables=dot_jac_statevariables)
         return
      end if
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! change of void ratio `\dot{e} = \left(1+e_i\right)\tr{(\mathbf{D})}`

      trT = Trace(cur_stress)                                        ! Trace of the stress tensor `\tr{(\mathbf{T})}`
      T_dev = Deviatoric_Part(cur_stress)                            ! Deviatoric part of the stress tensor `\mathbf{T}^\mathrm{dev}`
      T_dless = Dimensionless(cur_stress)                            ! Dimensionless stress tensor `\hat{\mathbf{T}}`
      T_dless_dev = Deviatoric_Part(T_dless)                         ! Deviatoric part of dimensionless stress tensor `\hat{\mathbf{T}}^\mathrm{dev}`

      ! Compute tensor `\mathcal{L}` and matrix `\mathbf{N}`
      factor_F = this%Get_Factor_F(stress_dless_dev=T_dless_dev)     ! (2.66) of Niemunis (2003): `F=\sqrt{\frac{1}{8}\tan^2{(\psi)} + \frac{2-\tan^2{(\psi)}}{2+\sqrt{2}\tan{(\psi)}\cos{(3\theta)}}} - \frac{1}{2\sqrt{2}} \tan{(\psi)}`
      factor_F2 = factor_F**2

      L_hat = (factor_F2 + this%derived_b2)*const_identity4d_sym &   ! (2.63) of Niemunis (2003) with increased shear stiffness from p.167:
            + this%derived_a2*Dyadic_Product22(T_dless, T_dless) &   ! `\hat{\mathcal{L}}_\mathrm{incr} = \left(F^2+b^2\right) \mathcal{I} + a^2\hat{\mathbf{T}}\otimes\hat{\mathbf{T}} - \frac{b^2}{3}\mathbf{I}\otimes\mathbf{I}`
            - this%derived_b2/3.0_dp*const_identity4d_tr

      L_hat_inv = 1.0_dp/factor_F2*(const_identity4d_sym &           ! (4.140) of Niemunis (2003): `\hat{\mathcal{L}}^{-1} = \frac{1}{F^2}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right)`
                - Dyadic_Product22(T_dless, T_dless)/(factor_F2/this%derived_a2 &
                + Double_Contraction22(T_dless, T_dless)))
      N_hat = factor_F*this%derived_a*(T_dless + T_dless_dev)
      B_mat = Double_Contraction42(L_hat_inv, N_hat)                 ! (4.141) of Niemunis (2003): `\mathbf{B} = \hat{\mathcal{L}}^{-1}:\hat{\mathbf{N}} = \frac{a}{F}\left(\mathcal{I}-\frac{\hat{\mathbf{T}}\otimes\hat{\mathbf{T}}}{\left(\frac{F}{a}\right)^2+\hat{\mathbf{T}}:\hat{\mathbf{T}}}\right):\left(\hat{\mathbf{T}}+\hat{\mathbf{T}}^\mathrm{dev}\right)`

      ! Alternatively used
      ! f_1 = factor_F2 + this%derived_a2*Double_Contraction22(T_dless, T_dless)
      ! f_2 = -this%derived_a2*Double_Contraction22(T_dless, T_dless+T_dless_dev)
      ! B_mat = f_1*(T_dless + T_dless_dev) + f_2 * T_dless

      B_dir = Nonzero_Division(val=B_mat, fac=Norm(B_mat))

      p_mean = -trT/3.0_dp
      q_mean = sqrt(1.5_dp*Double_Contraction22(T_dev, T_dev))

      crit_stress_ratio = factor_F * this%derived_pq_incl            ! p.126 of Niemunis (2003): `M(\mathbf{T}) = F(\mathbf{T}) \frac{6\sin{(\varphi_c)}}{3-\sin{(\varphi_c)}}`
      eta = q_mean/(crit_stress_ratio*p_mean)                        ! `\bar{\eta} = \frac{q}{Mp}`
      eta2 = eta**2                                                  ! Within (4.80) of Niemunis (2003): `\bar{\eta}^2 =\left(\frac{q}{Mp}\right)^2`

      pressure_preconsol = this%Get_Pressure_Preconsol( &            ! `p_e^+`
                           overcrit=eta2, trT=trT)
      state_OCR = pressure_equiv / pressure_preconsol                ! (4.78) of Niemunis (2003): `OCR = \frac{p_e}{p_e^+}`
      L_mat = -this%derived_beta_b*trT*L_hat                         ! (4.64) and (4.66) of Niemunis (2003): `\mathcal{L} = f_b\hat{\mathcal{L}} = -\beta_b\tr{(\mathbf{T})}\hat{\mathcal{L}}`

      D_vis = -this%param_D_r*B_dir &                                ! (4.74) of Niemunis (2003): `\mathbf{D}^\mathrm{vis} = -D_r \vec{\mathbf{B}}\left(\frac{1}{OCR}\right)^{(1/I_v)}`
            * (1.0_dp/state_OCR)**(1.0_dp/this%param_I_v)

      if (eta > 1.0_dp) then                                         ! State above critical state surface `||\mathbf{B}|| = 1`
         pressure_preconsol = p_mean*(1.0_dp + eta2)/2.0_dp &        ! (4.118) of Niemunis (2003): `p_e^{+\mathrm{new}} = p\left(1+\bar{\eta}^2\right)\frac{1+\beta_B}{2}`
                            * (1.0_dp + this%param_beta_b)
         dpressure_preconsol_dp = (1.0_dp - eta2)/2.0_dp &           ! Modification of (4.119) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial p} = \frac{p_e^{+\mathrm{new}}}{p}`
                                * (1.0_dp + this%param_beta_b)
         dpressure_preconsol_dq = eta/crit_stress_ratio &            ! Modification of (4.120) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial q} = \frac{\bar{\eta}(\beta_B + 1)}{M}`
                                * (1.0_dp + this%param_beta_b)
      else                                                           ! State below critical state surface `||\mathbf{B}|| = 1`
         pressure_preconsol = p_mean*(1.0_dp - this%param_beta_b &   ! (4.117) of Niemunis (2003): `p_e^{+\mathrm{new}} = \frac{p}{\beta_B - 1}\left[\beta_B\sqrt{1+\bar{\eta}^2(\beta_B^2 - 1)} - 1\right]`
                            * sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp)))/(1.0_dp - this%param_beta_b)
         dpressure_preconsol_dp = pressure_preconsol/p_mean &        ! (4.119) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial p} = \frac{p_e^{+\mathrm{new}}}{p} - \frac{\beta_B \bar{\eta}^2(\beta_B + 1)}{\sqrt{1 + \bar{\eta}^2(\beta_B^2-1)}}`
                                - this%param_beta_b*eta2 * (this%param_beta_b + 1.0_dp) &
                                / sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp))
         dpressure_preconsol_dq = this%param_beta_b*eta &            ! (4.120) of Niemunis (2003): `\frac{\partial p_e^{+\mathrm{new}}}{\partial q} = \frac{\beta_B \bar{\eta}(\beta_B + 1)}{M\sqrt{1 + \bar{\eta}^2(\beta_B^2-1)}}`
                                / crit_stress_ratio*(this%param_beta_b + 1.0_dp) &
                                / sqrt(1.0_dp + eta2*(this%param_beta_b**2 - 1.0_dp))
      end if

      ! Calculate `\frac{\partial p_e^{+\mathrm{new}}}{\partial \mathbf{T}}` as in (4.133) of Niemunis (2003)
      dpressure_preconsol_ds = -const_identity2d*dpressure_preconsol_dp/3.0_dp
      if (q_mean > 0.0001_dp) then
         dpressure_preconsol_ds = dpressure_preconsol_ds + T_dev*1.5_dp/q_mean*dpressure_preconsol_dq
      end if

      ! NOTE: Check dependence of dtth on ref_dt and compare it with other solutions. Although using dtth=0 would make
      !       the following procedure easier (L_mat and D_vis will stay), it might be critical to have a time dependence
      dtth = 0.5*ref_dt

      A_term = Dyadic_Product22(D_vis, dpressure_preconsol_ds)*dtth &! Compare to (4.132) of Niemunis (2003):
             * state_OCR/(this%param_I_v*pressure_equiv)             ! `2\mathcal{A} \approx \mathbf{D}^\mathrm{vis} \frac{OCR}{I_v \cdot p_e}\frac{\partial p_e^{+\mathrm{new}}}{\partial \mathbf{T}}`
      B_term = Dyadic_Product22(D_vis, const_identity2d)*dtth &      ! Compare to (4.134) of Niemunis (2003):
             * (1.0_dp + cur_voidratio) &                            ! `2\mathcal{B} = \frac{1+e^t}{I_v \lambda}\mathbf{D}_t^\mathrm{vis} \mathbf{1}`
             / (this%param_I_v*this%param_lambda)
      C_inv = const_identity4d_sym - B_term                          ! Inverse of (4.130) of Niemunis (2003): `\mathcal{C} = \left[\mathcal{M}^t-\mathcal{L}^t:\mathcal{B}\Delta t\right]^{-1}:\mathcal{L}^t`

      call Inverse_Tensor(tensor=C_inv, inv_tensor=C_mat, success=success)
      if (.not. success) then
         call Write_Error_And_Exit('Viscohypoplasticity: Inversion of C_inv failed')
      end if

      K_firstinv = Double_Contraction44(L_mat, A_term) &             ! Inverse of first part of (4.129) of Niemunis (2003): `\mathcal{I} + \mathcal{L}^t:\mathcal{A}\Delta t`
                 + const_identity4d_sym

      call Inverse_Tensor(tensor=K_firstinv, inv_tensor=K_first, success=success)
      if (.not. success) then
         call Write_Error_And_Exit('Viscohypoplasticity: Inversion of K_firstinv failed')
      end if

      L_mat = Double_Contraction44( &                                ! Similar to (4.129) of Niemunis (2003):
              Double_Contraction44(K_first, L_mat), C_inv)           ! `\mathcal{K} = \left[\mathcal{I} + \mathcal{L}^t:\mathcal{A}\Delta t\right]^{-1}:\left[\mathcal{M}^t - \mathcal{L}^t:\mathcal{B}\Delta t\right]`

      D_vis = Double_Contraction42(C_mat, D_vis)                     ! Apply transformation on `\mathbf{D}^\mathrm{vis}` as as below (4.131) of Niemunis (2003)

      N_mat = 0.0_dp
      call this%Intergranular_Strain(L_mat=L_mat, fdN_mat=N_mat, D_vis=D_vis, hypoplastic=.False., &
         igran_strain=igran_strain, R_max=this%param_R_max, m_T=this%param_m_T, m_R=this%param_m_R, &
         beta_R=this%param_beta_R, chi=this%param_chi, dt=ref_dt, dot_strain=dot_strain, &
         dot_stress=dot_stress, dot_igran_strain=dot_igran_strain, jacobian=dot_jac_stress)

      ! NOTE: Estimation of current stiffness for replacement model should be checked for correctness
      if ((this%calculate_jacobian) .and. (.not. setting_numerical_jacobian)) then
         cur_param_young = Tensor_Partialtrace(dot_jac_stress)/3.0_dp
      else
         cur_param_young = Tensor_Partialtrace(L_mat)/3.0_dp
      end if

      this%direct_variables(11) = cur_param_young

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_statevariables(2:10) = Mat_To_Vec9(mat=dot_igran_strain)
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Viscohypoplasticity_Ni03_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Ko15_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Ko15
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_c2, param_c3, param_c4, param_c5, param_e_c0, param_e_min
      real(dp) :: derived_c1

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval
      !
      class(Barodesy_Ko15), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      real(dp) :: K_c

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)

      ! --- Parameters
      this%param_phi_c = params(1)                                   ! Friction angle `\varphi_c` (in radians)
      this%param_c2    = params(2)                                   ! Parameter `c_2`
      this%param_c3    = params(3)                                   ! Parameter `c_3`
      this%param_c4    = params(4)                                   ! Parameter `c_4`
      this%param_c5    = params(5)                                   ! Parameter `c_5`
      this%param_e_c0  = params(6)                                   ! Critical void ratio `e_{c0}`
      this%param_e_min = params(7)                                   ! Minimal void ratio `e_\text{min}`

      ! --- Check valid parameter range
      if (firstcall) then
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('e_min', this%param_e_min, [setting_epsilon, 20.0_dp])
      end if

      K_c = (1.0_dp - sin(this%param_phi_c)) &                       ! `K_c = \frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}`
          / (1.0_dp + sin(this%param_phi_c))

      this%derived_c1 = sqrt(2.0_dp/3.0_dp)*log(K_c)                 ! (32) of Kolymbas (2015): `c_1 = \sqrt{\frac{2}{3}}\log{(K_c)}`
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: Nonzero_Division, Norm, Trace, Matrix_Exponential, Pack_States, Unpack_States
      !
      class(Barodesy_Ko15), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, R_mat, R_dir, D_dir, stress_dir
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: cur_voidratio, dot_voidratio, norm_D, norm_T, dilatancy, B_fak, e_c, &
                  h_fac, f_fac, g_fac

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)

      norm_D = Norm(dot_strain)
      D_dir = Nonzero_Division(val=dot_strain, fac=norm_D)
      dilatancy = Trace(D_dir)

      norm_T = Norm(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)

      R_mat = -Matrix_Exponential(mat=this%derived_c1*D_dir &        ! (9) of Kolymbas (2015): `\mathbf{R} = -e^{c_1\mathbf{D}^0 e^{(c_2 \delta)}}`
            * exp(this%param_c2*dilatancy))

      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))           ! `\mathbf{R}^0 = \frac{\mathbf{R}}{||\mathbf{R}||}`

      B_fak = (this%param_e_c0 - this%param_e_min) &                 ! (27) of Kolymbas (2015): `B = \frac{e_{c0} - e_\text{min}}{e_{c0} + 1} \left(\frac{c_4 + c_5\||\mathbf{T}||}{c_4}\right)^{-(1 + e_\text{min})/c_5}`
            / (this%param_e_c0 + 1.0_dp) &
            * ((this%param_c4 + this%param_c5*norm_T)/this%param_c4)**(-(1.0_dp + this%param_e_min)/this%param_c5)
      e_c = (this%param_e_min + B_fak)/(1.0_dp - B_fak)              ! (26) of Kolymbas (2015): `e_c = \frac{e_\text{min}+B}{1-B}`

      f_fac = dilatancy + this%param_c3*e_c                          ! (29) of Kolymbas (2015): `f = \delta + c_3 e_c`
      g_fac = -this%param_c3*cur_voidratio                           ! (30) of Kolymbas (2015): `g = -c_3 e`
      h_fac = -(this%param_c4 + this%param_c5*norm_T) &              ! (23) of Kolymbas (2015): `h = \frac{c_4 + c_5||\mathbf{T}||}{e - e_\text{min}}`
            / (cur_voidratio - this%param_e_min)

      dot_stress = h_fac*(f_fac*R_dir + g_fac*stress_dir)*norm_D     ! (12) of Kolymbas (2015): `\overset{\circ}{\mathbf{T}} = h\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Barodesy_Ko15_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Sc18_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Sc18
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_K_r, param_xi, param_kappa, param_e_c0, param_ce, &
                  param_c5, param_c6, param_c7
      real(dp) :: derived_c1, derived_c2, derived_c3, derived_c4, derived_c8, &
                  derived_zeta, derived_b_comp, derived_b_ext

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon
      use Math_Operations, only: Abort_If_Not_In_Interval, const_root2, const_root3
      !
      class(Barodesy_Sc18), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      ! Auxiliary terms `f_1`, `f_2`, `f_3` and `f_4`
      real(dp), parameter :: fak1 = log((1.0_dp + const_root2)/const_root2)
      real(dp), parameter :: fak2 = log((3.0_dp*const_root2 - const_root3)/(3.0_dp*const_root2))
      real(dp), parameter :: fak3 = log((1.0_dp + const_root2)/(2.0_dp + const_root2))
      real(dp), parameter :: fak4 = fak1*log((3.0_dp + 3.0_dp*const_root2) &
                                  / (3.0_dp + 3.0_dp*const_root2 - const_root3)) - fak3*fak2
      ! NOTE: Consider making cmin a class variable since it is also used in Calculate_Dot_State()
      real(dp), parameter :: cmin = -30.0_dp
      real(dp), parameter :: p_r = 1.0_dp                            ! Reference pressure `p_r` for parameters of barodey is 1 kPa
      real(dp) :: K_c, alpha_p, alpha_c, alpha_0

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)

      ! --- Parameters
      this%param_phi_c = params(1)                                   ! Friction angle `\varphi_c` (in radians)
      this%param_K_r   = params(2)                                   ! Reference stiffness `K_r` at 1 kPa
      this%param_xi    = params(3)                                   ! Exponent `\xi` in stiffness term
      this%param_kappa = params(4)                                   ! Stiffness releation `\kappa` between loading and unloading
      this%param_e_c0  = params(5)                                   ! Critical void ratio for zero stress `e_{c0}`
      this%param_ce    = params(6)                                   ! Ratio `c_e` of crit. and comparable void ratio at normal compression
      this%param_c5    = params(7)                                   ! Weighting factor `c_5` for pycnotropy
      this%param_c6    = params(8)                                   ! Control factor `c_6` for splitting `f` and `g`
      this%param_c7    = params(9)                                   ! Interpolation factor `c_7` of `b` between compression and extension

      ! --- Check valid parameter range
      if (firstcall) then
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
         call Abort_If_Not_In_Interval('K_r', this%param_K_r, [0.1_dp, 1.0e9_dp])
         call Abort_If_Not_In_Interval('xi', this%param_xi, [0.2_dp, 2.0_dp])
         call Abort_If_Not_In_Interval('kappa', this%param_kappa, [setting_epsilon, 1.0e9_dp])
         call Abort_If_Not_In_Interval('e_c0', this%param_e_c0, [setting_epsilon, 20.0_dp])
         call Abort_If_Not_In_Interval('c_e', this%param_ce, [1.0_dp, 5.0_dp])
         call Abort_If_Not_In_Interval('c_5', this%param_c5, [0.1_dp, 100.0_dp])
         call Abort_If_Not_In_Interval('c_6', this%param_c6, [0.1_dp, 100.0_dp])
         call Abort_If_Not_In_Interval('c_7', this%param_c7, [0.1_dp, 100.0_dp])
      end if

      ! --- Derived parameters (taken/adapted from F. Schranz (2018), p.166 and p.175)
      K_c = (1.0_dp - sin(this%param_phi_c)) &                       ! `K_c = \frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}`
          / (1.0_dp + sin(this%param_phi_c))

      alpha_p = log((9.0_dp - 11.0_dp*sin(this%param_phi_c)) &       ! `\alpha_p = \log{\left(\frac{9-11\sin{(\varphi_c)}}{9+13\sin{(\varphi_c)}}\right)}`
              / (9.0_dp + 13.0_dp*sin(this%param_phi_c)))*const_root3/2.0_dp
      alpha_c = const_root2/const_root3*log(K_c)                     ! `\alpha_c = \sqrt{\frac{2}{3}}\log{\left(\frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}\right)}`
      alpha_0 = log(1.0_dp - sin(this%param_phi_c))                  ! `\alpha_0 = \log{\left(1-\sin{(\varphi_c)}\right)}`

      this%derived_c1 = (log((alpha_p - cmin) / (alpha_c - cmin)) &  ! `c_1 = \frac{f_1}{f_4}\log{\left(\frac{\alpha_p+30}{\alpha_c+30}\right)} - \frac{f_2}{f_4}\log{\left(\frac{\alpha_0+30}{\alpha_c+30}\right)}`
                      * fak1 - log((alpha_0 - cmin)/(alpha_c - cmin)) * fak2)/fak4
      this%derived_c2 = (log((alpha_0 - cmin)/(alpha_c - cmin)) &    ! `c_2 = \frac{\log{\left(\frac{\alpha_0+30}{\alpha_c+30}\right)} - f_3 c_1}{f_1}`
                      - fak3*this%derived_c1)/fak1
      this%derived_c3 = (alpha_c - cmin) &                           ! `c_3 = \left(\alpha_c + 30\right)\frac{\left(1+\sqrt{2}\right)^{c_1}}{\sqrt{2}^{c_2}}`
                      * (1.0_dp + const_root2)**this%derived_c1 / const_root2**this%derived_c2
      this%derived_c4 = this%param_K_r &                             ! `c_4 = K_r\left(\sqrt{3} p_r\right)^{-\xi}`
                      * (const_root3*p_r)**(-this%param_xi)
      this%derived_zeta = -1.0_dp/K_c                                ! `\zeta = -\frac{1}{K_c} = -\frac{1+\sin{(\varphi_c)}}{1-\sin{(\varphi_c)}}`

      this%derived_b_comp = (this%param_c5 &                         ! (4.50) of Schranz (2018): `b_\mathrm{comp} = \frac{c_5\left(c_e^\zeta - 1\right) -3}{\sqrt{3}}`
                          * (this%param_ce**this%derived_zeta - 1.0_dp) - 3.0_dp)/const_root3
      this%derived_b_ext = (this%param_c5 &                          ! (4.56) of Schranz (2018): `b_\mathrm{ext} = \frac{c_5\left(1 - c_e^\zeta\right) -3\kappa}{\sqrt{3}}`
                         * (1.0_dp - this%param_ce**this%derived_zeta) - 3.0_dp*this%param_kappa)/const_root3

      this%derived_c8 = -this%param_c6 &                             ! (4.65) of Schranz (2018): `c_8 = \frac{-c_6}{\sqrt{3} b_\mathrm{ext}}`
                      / (const_root3*this%derived_b_ext)
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: const_root3, Nonzero_Division, Norm, Trace, Matrix_Exponential, &
                                 Pack_States, Unpack_States
      !
      class(Barodesy_Sc18), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, R_mat, R_dir, D_dir, stress_dir
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: cur_voidratio, dot_voidratio, alpha, norm_D, norm_T, dilatancy, trT, e_c, &
                  h_fac, b_interp, f_fac, g_fac

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)

      norm_D = Norm(dot_strain)
      D_dir = Nonzero_Division(val=dot_strain, fac=norm_D)
      dilatancy = Trace(D_dir)                                       ! (1.7) of Schranz (2018): `\delta = \frac{\tr{(\mathbf{D})}}{||\mathbf{D}||}`

      norm_T = Norm(cur_stress)
      trT = Trace(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)

      alpha = -30.0_dp + this%derived_c3 &                           ! `\alpha = -30 + c_3\frac{|\delta - \sqrt{2}|^{c_2}}{\left(1+|\delta-\sqrt{2}|\right)^{c_1}}`
            * abs(dilatancy - sqrt(2.0_dp))**this%derived_c2/(1.0_dp + abs(dilatancy - sqrt(2.0_dp)))**this%derived_c1
      R_mat = -Matrix_Exponential(mat=D_dir*alpha)                   ! `\mathbf{R} = -e^{\alpha\mathbf{D}^0}`

      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))           ! `\mathbf{R}^0 = \frac{\mathbf{R}}{||\mathbf{R}||}`

      e_c = (1.0_dp + this%param_e_c0) &                             ! Modifed version of Appendix C of Schranz (2018):
          * exp(-(-trT/3.0_dp)**(1.0_dp - this%param_xi) &           ! `e_c = \left(1+e_{c0}\right)e^{\left(-\frac{p^{1-\xi}}{K_r\left(1-\xi\right)}\right)} - 1`
          / (this%param_K_r*(1.0_dp - this%param_xi))) - 1.0_dp

      h_fac = this%derived_c4*norm_T**this%param_xi                  ! (4.32) of Schranz (2018): `h = c_4 ||\mathbf{T}||^\xi`

      b_interp = (this%derived_b_ext - this%derived_b_comp) &        ! (4.57) of Schranz (2018): `b = \left(b_\mathrm{ext}-b_\mathrm{comp}\right)\left(\frac{\delta+\sqrt{3}}{2\sqrt{3}}\right)^{c_7} + b_\mathrm{comp}`
               * ((dilatancy + const_root3)/(2.0_dp*const_root3))**this%param_c7 + this%derived_b_comp

      f_fac = this%derived_c8*b_interp*dilatancy + this%param_c6     ! (4.59) of Schranz (2018): `f = c_8 b \delta + c_6`
      g_fac = (1.0_dp - this%derived_c8)*b_interp*dilatancy &        ! (4.59) of Schranz (2018): `g = (1-c_8) b \delta + c_5\left[\left(\frac{1+e}{1+e_c}\right)^\zeta - 1\right] -c_6`
            + this%param_c5*(((1.0_dp + cur_voidratio)/(1.0_dp + e_c))**this%derived_zeta - 1.0_dp) - this%param_c6

      dot_stress = h_fac*(f_fac*R_dir + g_fac*stress_dir)*norm_D     ! (3.14) of Schranz (2018): `\overset{\circ}{\mathbf{T}} = h\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)     ! As (4.33) of Schranz (2018)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Barodesy_Sc18_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Barodesy_Ko21_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Barodesy_Ko21
   ! --------------------------------------------------------------- !
      private
      real(dp) :: param_phi_c, param_c2, param_c3, param_c4, param_c5, param_k1, param_k2, &
                  param_kappa1, param_kappa2, param_e_c0
      real(dp) :: derived_c1

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: Abort_If_Not_In_Interval
      !
      class(Barodesy_Ko21), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      real(dp) :: K_c

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)
      this%direct_variables_mask(3:12) = .True.

      ! --- Parameters
      this%param_phi_c  = params(1)                                  ! Friction angle `\varphi_c` (in radians)
      this%param_c2     = params(2)                                  ! Parameter `c_2`
      this%param_c3     = params(3)                                  ! Parameter `c_3`
      this%param_c4     = params(4)                                  ! Parameter `c_4`
      this%param_c5     = params(5)                                  ! Parameter `c_5`
      this%param_kappa1 = params(6)                                  ! Parameter `\kappa_1`
      this%param_kappa2 = params(7)                                  ! Parameter `\kappa_2`
      this%param_k1     = params(8)                                  ! Parameter `k_1`
      this%param_k2     = params(9)                                  ! Parameter `k_2`
      this%param_e_c0   = params(10)                                 ! Critical void ratio `e_{c0}`

      ! --- Check valid parameter range
      if (firstcall) then
         call Abort_If_Not_In_Interval('phi_c', this%param_phi_c, [0.1745_dp, 1.0472_dp])
      end if

      K_c = (1.0_dp - sin(this%param_phi_c)) &                       ! `K_c = \frac{1-\sin{(\varphi_c)}}{1+\sin{(\varphi_c)}}`
          / (1.0_dp + sin(this%param_phi_c))
      this%derived_c1 = sqrt(2.0_dp/3.0_dp)*log(K_c)                 ! (16.26) of Kolymbas (2022): `c_1 = \sqrt{\frac{2}{3}}\log{(K_c)}`
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: const_identity2d, Nonzero_Division, Norm, Trace, Deviatoric_Part, &
                                 Matrix_Exponential, Pack_States, Unpack_States, Vec9_To_Mat, Mat_To_Vec9
      !
      class(Barodesy_Ko21), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress, stress_dir, R_mat, R_dir
      real(dp), dimension(6) :: D_dir, D_dirold, D_dirdev
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      real(dp) :: cur_voidratio, dot_voidratio, p_mean, norm_D, norm_T, delta, e_c, lambda_D, &
                  dot_e_c1, dot_e_c2, eps1, eps2, h_fac, f_fac, g_fac

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)
      cur_voidratio = statevariables(1)
      ! Currently e_c is saves in two fields due to different integration of both components
      e_c = statevariables(2) + statevariables(3)
      D_dirold = Vec9_To_Mat(vec9=statevariables(4:12))

      norm_D = Norm(dot_strain)
      D_dir = Nonzero_Division(val=dot_strain, fac=norm_D)
      delta = Trace(D_dir)

      norm_T = Norm(cur_stress)
      stress_dir = Nonzero_Division(val=cur_stress, fac=norm_T)
      p_mean = -Trace(cur_stress)/3.0_dp

      eps1 = this%param_k1 - this%param_k2*log(p_mean)               ! (1) of Kolymbas (2021): `\epsilon_1(p) = k_1 - k_2 \log{(p)}`
      eps2 = eps1*(1.0_dp + this%param_kappa2*delta)                 ! (2) of Kolymbas (2021): `\epsilon_2(p) = \epsilon_1(p)\cdot \left(1 + \kappa_2 \delta\right)`

      D_dirdev = Deviatoric_Part(D_dir)
      R_mat = -Matrix_Exponential(mat=D_dirdev * this%derived_c1) &  ! (11) of Kolymbas (2021): `\mathbf{R}(\mathbf{D}) = -\exp{\left(c_1 \hat{\mathbf{D}}^0\right)} + c_2 \delta \mathbf{1}`
            + this%param_c2*delta*const_identity2d
      R_dir = Nonzero_Division(val=R_mat, fac=Norm(R_mat))

      h_fac = this%param_c4 * norm_T**this%param_c5                  ! (12) of Kolymbas (2021): `h = c_4 ||\mathbf{T}||^{c_5}`
      f_fac = e_c + this%param_c3*delta                              ! `f = e_c + c_3 \delta`
      g_fac = -cur_voidratio + this%param_c3*delta                   ! `g = -e + c_3 \delta`

      dot_stress = h_fac * (f_fac*R_dir + g_fac*stress_dir) * norm_D ! (10) of Kolymbas (2021): `\overset{\circ}{\mathbf{T}} = h\cdot\left(f\mathbf{R}^0 + g\mathbf{T}^0\right)||\mathbf{D}||`
      dot_voidratio = (1.0_dp + cur_voidratio)*Trace(dot_strain)

      ! The integration of dot_e_c relies on two different parts as hinted in Kolymbas (2021)
      dot_e_c1 = this%param_kappa1*(eps2 - e_c)                      ! First part of (3) of Kolymbas (2021): `\dot{e}_{c,1} = \kappa_1 \cdot \left(\epsilon_2(p) - e_c\right)\cdot \dot{\varepsilon}`
      !                                                              ! which will be integrated normally
      lambda_D = Norm(D_dir - D_dirold)
      dot_e_c2 = lambda_D * 0.5_dp * (this%param_e_c0 - e_c)         ! Second part of (3) of Kolymbas (2021): `\dot{e}_{c,2} =  \left(e_{c0} - e_c\right) \cdot \dot{\lambda_D} / 2`
      !                                                              ! which has to be added to the integrated `\dot{e}_{c,1}`

      this%direct_variables(3) = dot_e_c2
      this%direct_variables(4:12) = Mat_To_Vec9(mat=D_dir)

      ! --- Packing all dot_states in a long vector
      dot_statevariables(1) = dot_voidratio
      dot_statevariables(2) = dot_e_c1
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Barodesy_Ko21_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Test_DGL_Class
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Constitutive_Model) :: Test_DGL
   ! --------------------------------------------------------------- !
      private
      integer :: param_select

      contains

      procedure :: Initialize
      procedure :: Calculate_Dot_State
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Test_DGL), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !

      ! --- Initialize base class variables
      call this%Base_Initialization(calculate_jacobian=calculate_jacobian, provide_jacobian=.False.)

      ! --- Parameters
      this%param_select = int(params(1))                             ! Select, which test DGL should be used
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   function Calculate_Dot_State(this, ref_dt, cur_time, cur_state, dot_strain) result(dot_state)
   ! --------------------------------------------------------------- !
      use Math_Operations, only: Pack_States, Unpack_States, Ref_Index, Set_Element_In_Tensor, Double_Contraction42
      !
      class(Test_DGL), intent(inout) :: this
      real(dp), intent(in) :: ref_dt
      real(dp), intent(in) :: cur_time
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), dimension(6), intent(in) :: dot_strain
      real(dp), dimension(setting_max_internal_states) :: dot_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: cur_stress, dot_stress
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables
      real(dp), dimension(6, 6) :: jac_stress, dot_jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables, dot_jac_statevariables
      !
      real(dp) :: assign_val
      integer :: ref_idx, ref_jdx, idx, jdx, kdx, ldx
      real(dp) :: mod_idx, mod_jdx
      logical :: symmetric

      dot_state = 0.0_dp
      dot_jac_stress = 0.0_dp
      dot_statevariables = 0.0_dp
      dot_jac_statevariables = 0.0_dp

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=cur_state, stress=cur_stress, jac_stress=jac_stress, &
         statevariables=statevariables, jac_statevariables=jac_statevariables)

      ! --- Calculate dot_stress (response of test DGL) based on the selected test case
      if (this%param_select == 1) then
         dot_stress = 10.0_dp * cur_stress
         ! Solution for this DGL is:   stress_0 * exp(10*t)
      else if (this%param_select == 2) then
         dot_stress = 0.8_dp*(exp(2.5_dp)-exp(0.5_dp))*exp(-0.25_dp*cur_stress)
         ! Solution for this DGL is:   4 * log(0.2*t*(exp(2.5)-exp(0.5)) + exp(0.5))
      else if (this%param_select == 3) then
         dot_stress = log(9.0_dp)/5.0_dp*(cur_stress - 1.0_dp)
         ! Solution for this DGL is:   9**(0.2*t) + 1
      else
         symmetric = .True.
         do ref_idx = 1, 6
            call Ref_Index(ref_idx=ref_idx, idx=idx, jdx=jdx)
            do ref_jdx = 1, 6
               call Ref_Index(ref_idx=ref_jdx, idx=kdx, jdx=ldx)
               assign_val = 10.0_dp*(ref_idx-1)+ref_jdx
               if (6 == 6) then
                  mod_idx = real(ref_idx, dp)
                  if (ref_idx > 3) then
                     mod_idx = 2 + (mod_idx - 3)*2
                  end if
                  mod_jdx = real(ref_jdx, dp)
                  if (ref_jdx > 3) then
                     mod_jdx = 2 + (mod_jdx - 3)*2
                  end if
                  assign_val = 10.0_dp*(mod_idx-1)+mod_jdx
               else if (symmetric) then
                  mod_idx = real(ref_idx, dp)
                  if (ref_idx > 3) then
                     mod_idx = 2 + floor((mod_idx - 2)/2)*2
                  end if
                  mod_jdx = real(ref_jdx, dp)
                  if (ref_jdx > 3) then
                     mod_jdx = 2 + floor((mod_jdx - 2)/2)*2
                  end if
                  assign_val = 10.0_dp*(mod_idx-1)+mod_jdx
               end if
               call Set_Element_In_Tensor(tens=dot_jac_stress, idx=idx, jdx=jdx, kdx=kdx, ldx=ldx, val=assign_val)
            end do
            dot_stress = Double_Contraction42(dot_jac_stress, dot_strain)
         end do
      end if

      ! --- Packing all dot_states in a long vector
      dot_state = Pack_States(stress=dot_stress, jac_stress=dot_jac_stress, statevariables=dot_statevariables, &
         jac_statevariables=dot_jac_statevariables)
   end function Calculate_Dot_State
end module Test_DGL_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Constitutive_Model_Class
   use General_Settings, only: dp, Write_Error_And_Exit, setting_len_id, setting_id_elasticity, setting_id_hypo_wu92, &
                               setting_id_hypo_vw96, setting_id_viscohypo_ni03, setting_id_barodesy_ko15, &
                               setting_id_barodesy_sc18, setting_id_barodesy_ko21, setting_id_test_dgl
   use Constitutive_Model_Baseclass
   implicit none

   private
   public :: Select_Constitutive_Model, Constitutive_Model


   contains


   ! --------------------------------------------------------------- !
   subroutine Select_Constitutive_Model(new_constitutive_model, identifier, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use Elasticity_Class
      use Hypoplasticity_Wu92_Class
      use Hypoplasticity_VW96_Class
      use Viscohypoplasticity_Ni03_Class
      use Barodesy_Ko15_Class
      use Barodesy_Sc18_Class
      use Barodesy_Ko21_Class
      use Test_DGL_Class
      !
      class(Constitutive_Model), allocatable, intent(out) :: new_constitutive_model
      character(len=setting_len_id), intent(in) :: identifier
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      if (identifier == setting_id_elasticity) then
         allocate(Elasticity::new_constitutive_model)
      else if (identifier == setting_id_hypo_wu92) then
         allocate(Hypoplasticity_Wu92::new_constitutive_model)
      else if (identifier == setting_id_hypo_vw96) then
         allocate(Hypoplasticity_VW96::new_constitutive_model)
      else if (identifier == setting_id_viscohypo_ni03) then
         allocate(Viscohypoplasticity_Ni03::new_constitutive_model)
      else if (identifier == setting_id_barodesy_ko15) then
         allocate(Barodesy_Ko15::new_constitutive_model)
      else if (identifier == setting_id_barodesy_sc18) then
         allocate(Barodesy_Sc18::new_constitutive_model)
      else if (identifier == setting_id_barodesy_ko21) then
         allocate(Barodesy_Ko21::new_constitutive_model)
      else if (identifier == setting_id_test_dgl) then
         allocate(Test_DGL::new_constitutive_model)
      else
         call Write_Error_And_Exit('Select_Constitutive_Model: Identifier >' // identifier // '< unknown')
      end if
      call new_constitutive_model%Initialize(params=params, calculate_jacobian=calculate_jacobian, &
         firstcall=firstcall)
   end subroutine Select_Constitutive_Model
end module Constitutive_Model_Class


! ==================================================================================================================== !
module Solver_Baseclass
   use General_Settings, only: dp, setting_num_statevariables, setting_max_internal_states
   use Constitutive_Model_Baseclass, only: Constitutive_Model

   implicit none

   private
   public :: Solver

   ! --------------------------------------------------------------- !
   type, abstract :: Solver
   ! --------------------------------------------------------------- !
      character(len=24) :: name
      logical :: current_error_assigned, next_dot_state_assigned, stepsize_fixed, last_integ_success
      integer :: last_num_steps_accepted, last_num_steps_rejected
      real(dp), dimension(setting_max_internal_states) :: current_error, next_dot_state
      logical, dimension(setting_max_internal_states) :: mask_ignored_indices
      real(dp) :: exp_stepgrow, max_stepgrow, exp_stepshrink, min_stepshrink

      contains

      procedure(integration_algorithm_interface), deferred :: Integration_Algorithm

      procedure :: Initialize
      procedure :: Integrate
      procedure :: Set_State_Mask
      procedure :: Set_Error_Estimate
      procedure :: Get_Statistics
   end type

   ! --------------------------------------------------------------- !
   abstract interface
   ! --------------------------------------------------------------- !
      function integration_algorithm_interface(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) &
         result(new_state)
         import
         class(Solver), intent(inout) :: this
         class(Constitutive_Model), intent(inout) :: integ_obj
         real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
         real(dp), intent(in) :: cur_time, dt
         real(dp), dimension(setting_max_internal_states) :: new_state
      end function
   end interface


   contains


   ! --------------------------------------------------------------- !
   subroutine Initialize(this, name, exp_stepgrow, max_stepgrow, exp_stepshrink, min_stepshrink, stepsize_fixed)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), optional :: exp_stepgrow
      real(dp), intent(in), optional :: max_stepgrow
      real(dp), intent(in), optional :: exp_stepshrink
      real(dp), intent(in), optional :: min_stepshrink
      logical, intent(in), optional :: stepsize_fixed
      ! ------------------------------------------------------------ !
      this%name = name

      this%current_error_assigned = .False.
      this%current_error = 0.0_dp

      this%next_dot_state_assigned = .False.
      this%next_dot_state = 0.0_dp

      this%mask_ignored_indices = .False.

      this%last_integ_success = .False.
      this%last_num_steps_accepted = 0
      this%last_num_steps_rejected = 0

      if (present(exp_stepgrow)) then
         this%exp_stepgrow = exp_stepgrow
      else
         this%exp_stepgrow = -0.2_dp
      end if

      if (present(exp_stepshrink)) then
         this%exp_stepshrink = exp_stepshrink
      else
         this%exp_stepshrink = -0.25_dp
      end if

      if (present(max_stepgrow)) then
         this%max_stepgrow = max_stepgrow
      else
         this%max_stepgrow = 5.0_dp
      end if

      if (present(min_stepshrink)) then
         this%min_stepshrink = min_stepshrink
      else
         this%min_stepshrink = 0.1_dp
      end if

      if (present(stepsize_fixed)) then
         this%stepsize_fixed = stepsize_fixed
      else
         this%stepsize_fixed = .False.
      end if
   end subroutine Initialize


   ! --------------------------------------------------------------- !
   subroutine Integrate(this, integ_obj, cur_state, cur_time, dt, new_state, mod_dt, successful_integration)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_epsilon, Write_Warning, &
                                  setting_stepsize_scaling_safety, setting_max_rel_error, &
                                  setting_max_integration_steps, setting_max_integration_refinements
      use Debug, only: Formatval
      !
      class(Solver), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states), intent(out) :: new_state
      real(dp), intent(inout) :: mod_dt
      logical, intent(out) :: successful_integration
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k1, inp_state, cum_error, state_scaled_err
      real(dp), dimension(setting_num_statevariables) :: tmp_direct_variables, direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask
      real(dp) :: steptime, stepsize, temp_error, temp_stepsize, next_stepsize
      integer :: idx_integloop, idx_convloop

      successful_integration = .False.
      this%last_integ_success = .False.
      this%last_num_steps_accepted = 0
      this%last_num_steps_rejected = 0
      cum_error = 0.0_dp
      inp_state = cur_state

      direct_variables = 0.0_dp

      idx_integloop = 0
      idx_convloop = 0

      steptime = 0.0_dp
      stepsize = mod_dt
      if (stepsize > dt) then
         ! Initial mod_dt should not be larger than overall dt
         stepsize = dt
      end if

      if ((this%stepsize_fixed) .and. (stepsize < dt/setting_max_integration_steps)) then
         call Write_Warning('Integrate: Fixed step integrator needs more steps (' // &
            Formatval('', int(dt/stepsize)) // ') than allowed ' &
            // 'in setting_max_integration_steps (' // Formatval('', setting_max_integration_steps) // ')')
         new_state = cur_state
         return
      end if

      integration_loop: &
      do
         idx_integloop = idx_integloop + 1

         if (this%next_dot_state_assigned) then
            dot_state_k1 = this%next_dot_state
            this%next_dot_state_assigned = .False.
            this%next_dot_state = 0.0_dp
         else
            dot_state_k1 = integ_obj%Get_Dot_State( &                ! `k_1 = f(t_n, y_n)`
               ref_dt=stepsize, cur_time=cur_time+steptime, cur_state=inp_state)
         end if
         call integ_obj%Get_Direct_Variables(direct_variables=tmp_direct_variables)

         ! Scaled error estimation from Press et al. (1997): Numerical Recipes in Fortran, p. 711/715
         state_scaled_err = abs(inp_state) + abs(stepsize*dot_state_k1) + setting_epsilon
         next_stepsize = stepsize

         idx_convloop = 0                                            ! Reset before entering convergence_loop
         convergence_loop: &
         do
            idx_convloop = idx_convloop + 1
            stepsize = next_stepsize

            new_state = this%Integration_Algorithm(integ_obj=integ_obj, cur_state=inp_state, &
               dot_state_k1=dot_state_k1, cur_time=cur_time + steptime, dt=stepsize)
            ! Check if stepsize is not fixed, error information is present and error is too big to reject step
            if (this%stepsize_fixed) then
               if (this%current_error_assigned) then
                  cum_error = cum_error + this%current_error
               end if

               this%last_num_steps_accepted = this%last_num_steps_accepted + 1
               exit convergence_loop
            end if

            if (.not. this%current_error_assigned) then
               this%last_num_steps_accepted = this%last_num_steps_accepted + 1
               exit convergence_loop
            end if

            temp_error = maxval(abs(this%current_error)/state_scaled_err)/setting_max_rel_error
            if (temp_error <= 1.0_dp) then
               ! Step is accepted and scaling is performed on next step size
               cum_error = cum_error + this%current_error
               this%last_num_steps_accepted = this%last_num_steps_accepted + 1

               next_stepsize = setting_stepsize_scaling_safety*stepsize*(temp_error**this%exp_stepgrow)
               next_stepsize = min(next_stepsize, this%max_stepgrow*stepsize)
               exit convergence_loop
            else
               ! Step is rejected. Step size for this step will be scaled
               this%last_num_steps_rejected = this%last_num_steps_rejected + 1

               temp_stepsize = setting_stepsize_scaling_safety*stepsize*(temp_error**this%exp_stepshrink)
               next_stepsize = max(temp_stepsize, this%min_stepshrink*stepsize)
            end if

            if (idx_convloop >= setting_max_integration_refinements) then
               call Write_Warning('Integrate: Too many iterations in convergence_loop')
               return
            end if

            if (next_stepsize <= setting_epsilon) then
               call Write_Warning('Integrate: Stepsize got reduced below minimum value')
               return
            end if
         end do convergence_loop

         direct_variables = direct_variables + tmp_direct_variables*stepsize/dt

         inp_state = new_state

         steptime = steptime + stepsize
         stepsize = next_stepsize
         mod_dt = stepsize

         if (steptime >= dt - setting_epsilon) then
            exit integration_loop
         end if

         if (steptime + stepsize > dt) then
            stepsize = dt - steptime
         end if

         if (idx_integloop >= setting_max_integration_steps) then
            call Write_Warning('Integrate: Too many iterations in integration_loop')
            return
         end if
      end do integration_loop

      if (this%current_error_assigned) then
         this%current_error = cum_error
      end if
      successful_integration = .True.
      this%last_integ_success = .True.

      ! Assign direct variables to their positions in new_state
      direct_variables_mask = integ_obj%Get_Direct_Variables_Mask()
      where (direct_variables_mask)
         new_state((6+1)*6+1:(6+1)*6+setting_num_statevariables) = direct_variables
      end where
   end subroutine Integrate


   ! --------------------------------------------------------------- !
   subroutine Set_State_Mask(this, statemask)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      logical, dimension(setting_num_statevariables), intent(in) :: statemask
      ! ------------------------------------------------------------ !
      this%mask_ignored_indices((6+1)*6+1:(6+1)*6+setting_num_statevariables) = statemask
   end subroutine Set_State_Mask


   ! --------------------------------------------------------------- !
   subroutine Set_Error_Estimate(this, estimation)
   ! --------------------------------------------------------------- !
      class(Solver), intent(inout) :: this
      real(dp), dimension(setting_max_internal_states), intent(in) :: estimation
      ! ------------------------------------------------------------ !
      ! Mask all values (set to zero) which should not contribute to error estimation
      where (this%mask_ignored_indices)
         this%current_error = 0.0_dp
      elsewhere
         this%current_error = estimation
      end where
      this%current_error_assigned = .True.
   end subroutine Set_Error_Estimate


   ! --------------------------------------------------------------- !
   function Get_Statistics(this) result(statistics)
   ! --------------------------------------------------------------- !
      class(Solver), intent(in) :: this
      integer, dimension(3) :: statistics
      ! ------------------------------------------------------------ !
      integer :: integ_success

      integ_success = 0
      if (this%last_integ_success) then
         integ_success = 1
      end if
      statistics = [integ_success, this%last_num_steps_accepted, this%last_num_steps_rejected]
   end function Get_Statistics
end module Solver_Baseclass


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Euler_Explicit_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: Euler_Explicit
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(Euler_Explicit), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      new_state = cur_state + dt*dot_state_k1                        ! `y_{n+1} = y_n + h f(t_n, y_n)`
   end function Integration_Algorithm
end module Euler_Explicit_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Euler_Richardson_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: Euler_Richardson
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(Euler_Richardson), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.5h, y_n+0.5h k_1)`
         cur_time=cur_time + 0.5_dp*dt, cur_state=cur_state + dt*0.5_dp*dot_state_k1)
      new_state = cur_state + dt*dot_state_k2                        ! `y_{n+1} = y_n + h k_2`

      extra_state = cur_state + dt*(0.5_dp*dot_state_k1 + 0.5_dp*dot_state_k2)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module Euler_Richardson_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK23_Simpson_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK23_Simpson
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK23_Simpson), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n + 0.5h, y_n + 0.5h k_1)`
         cur_time=cur_time + 0.5_dp*dt, cur_state=cur_state + dt*0.5_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f(t_n + h, y_n - h k_1 + 2h k_2)`
         cur_time=cur_time + dt, cur_state=cur_state + dt*(-dot_state_k1 + 2.0_dp*dot_state_k2))

      new_state = cur_state + dt*(dot_state_k1/6.0_dp &              ! `y_{n+1} = y_n + \frac{1}{6}h k_1 + \frac{2}{3}h k_2 + \frac{1}{6}h k_3`
                + 2.0_dp/3.0_dp*dot_state_k2 + dot_state_k3/6.0_dp)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k3

      extra_state = cur_state + dt*dot_state_k2                      ! `y_{p,n+1} = y_n + h k_2`
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module RK23_Simpson_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK23_Bogacki_Shampine_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK23_Bogacki_Shampine
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK23_Bogacki_Shampine), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.5h, y_n+0.5h k_1)`
         cur_time=cur_time + 0.5_dp*dt, cur_state=cur_state + dt*0.5_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f(t_n+0.75h, y_n+0.75h k_2)`
         cur_time=cur_time + 0.75_dp*dt, cur_state=cur_state + dt*0.75_dp*dot_state_k2)

      new_state = cur_state + dt*(2.0_dp/9.0_dp*dot_state_k1 &       ! `y_{n+1} = y_n + \frac{2}{9}h k_1 + \frac{1}{3}h k_2 + \frac{4}{9}h k_3`
                + 1.0_dp/3.0_dp*dot_state_k2 + 4.0_dp/9.0_dp*dot_state_k3)

      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f(t_n+h, y_{n+1})`
         cur_time=cur_time + dt, cur_state=new_state)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k4

      extra_state = cur_state + dt*(7.0_dp/24.0_dp*dot_state_k1 &    ! `y_{p,n+1} = y_n + \frac{7}{24}h k_1 + \frac{1}{4}h k_2 + \frac{1}{3}h k_3 + \frac{1}{8}h k_4`
                  + 1.0_dp/4.0_dp*dot_state_k2 + 1.0_dp/3.0_dp*dot_state_k3 + 1.0_dp/8.0_dp*dot_state_k4)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module RK23_Bogacki_Shampine_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK45_Cash_Karp_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK45_Cash_Karp
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Cash_Karp), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.2h, y_n+0.2h k_1)`
         cur_time=cur_time + 0.2_dp*dt, &
         cur_state=cur_state + dt*0.2_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{10}h, y_n+\frac{3}{40}h k_1 + \frac{9}{40}h k_2\right)`
         cur_time=cur_time + 0.3_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/40.0_dp*dot_state_k1 + 9.0_dp/40.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{3}{5}h, y_n+\frac{3}{10}h k_1 - \frac{9}{10}h k_2 + \frac{6}{5}h k_3\right)`
         cur_time=cur_time + 0.6_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/10.0_dp*dot_state_k1 - 9.0_dp/10.0_dp*dot_state_k2 &
            + 6.0_dp/5.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+h, y_n-\frac{11}{54}h k_1 + \frac{5}{2}h k_2 - \frac{70}{27}h k_3 + \frac{35}{27}h k_4\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(-11.0_dp/54.0_dp*dot_state_k1 + 5.0_dp/2.0_dp*dot_state_k2 &
            - 70.0_dp/27.0_dp*dot_state_k3 + 35.0_dp/27.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+\frac{7}{8}h, y_n+\frac{1631}{55296}h k_1 + \frac{175}{512}h k_2 + \frac{575}{13824}h k_3 + \frac{44275}{110592}h k_4 + \frac{253}{4096}h k_5\right)`
         cur_time=cur_time + 7.0_dp/8.0_dp*dt, &
         cur_state=cur_state + dt*(1631.0_dp/55296.0_dp*dot_state_k1 &
            + 175.0_dp/512.0_dp*dot_state_k2 + 575.0_dp/13824.0_dp*dot_state_k3 &
            + 44275.0_dp/110592.0_dp*dot_state_k4 + 253.0_dp/4096.0_dp*dot_state_k5))

      new_state = cur_state + dt*(37.0_dp/378.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{37}{378}h k_1 + \frac{250}{621}h k_3 + \frac{125}{594}h k_4 - \frac{512}{1771}h k_6`
                + 250.0_dp/621.0_dp*dot_state_k3 + 125.0_dp/594.0_dp*dot_state_k4 &
                + 512.0_dp/1771.0_dp*dot_state_k6)

      extra_state = cur_state &                                      ! `y_{p,n+1} =  y_n + \frac{2825}{27648}h k_1 + \frac{18575}{48384}h k_3 + \frac{13525}{55296}h k_4 + \frac{277}{14336}h k_5 + \frac{1}{4}h k_6`
                  + dt*(2825.0_dp/27648.0_dp*dot_state_k1 + 18575.0_dp/48384.0_dp*dot_state_k3 &
                  + 13525.0_dp/55296.0_dp*dot_state_k4 + 277.0_dp/14336.0_dp*dot_state_k5 + 0.25_dp*dot_state_k6)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module RK45_Cash_Karp_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK45_Dormand_Prince_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK45_Dormand_Prince
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Dormand_Prince), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, dot_state_k7, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.2h, y_n+0.2h k_1)`
         cur_time=cur_time + 0.2_dp*dt, &
         cur_state=cur_state + dt*0.2_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{10}h, y_n+\frac{3}{40}h k_1 + \frac{9}{40}h k_2\right)`
         cur_time=cur_time + 0.3_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/40.0_dp*dot_state_k1 + 9.0_dp/40.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{4}{5}h, y_n+\frac{44}{45}h k_1 - \frac{56}{15}h k_2 + \frac{32}{9}h k_3\right)`
         cur_time=cur_time + 0.8_dp*dt, &
         cur_state=cur_state + dt*(44.0_dp/45.0_dp*dot_state_k1 - 56.0_dp/15.0_dp*dot_state_k2 &
            + 32.0_dp/9.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+\frac{8}{9}h, y_n+\frac{19372}{6561}h k_1 - \frac{25360}{2187}h k_2 + \frac{64448}{6561}h k_3 - \frac{212}{729}h k_4\right)`
         cur_time=cur_time + 8.0_dp/9.0_dp*dt, &
         cur_state=cur_state + dt*(19372.0_dp/6561.0_dp*dot_state_k1 &
            - 25360.0_dp/2187.0_dp*dot_state_k2 + 64448.0_dp/6561.0_dp*dot_state_k3 &
            - 212.0_dp/729.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+h, y_n+\frac{9017}{3168}h k_1 - \frac{355}{33}h k_2 + \frac{46732}{5247}h k_3 + \frac{49}{176}h k_4 - \frac{5103}{18656}h k_5\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(9017.0_dp/3168.0_dp*dot_state_k1 - 355.0_dp/33.0_dp*dot_state_k2 &
            + 46732.0_dp/5247.0_dp*dot_state_k3 + 49.0_dp/176.0_dp*dot_state_k4 &
            - 5103.0_dp/18656.0_dp*dot_state_k5))

      new_state = cur_state + dt*(35.0_dp/384.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{35}{384}h k_1 + \frac{500}{1113}h k_3 + \frac{125}{192}h k_4 - \frac{2187}{6784}h k_5 + \frac{11}{84}h k_6`
                + 500.0_dp/1113.0_dp*dot_state_k3 + 125.0_dp/192.0_dp*dot_state_k4 &
                - 2187.0_dp/6784.0_dp*dot_state_k5 + 11.0_dp/84.0_dp*dot_state_k6)

      dot_state_k7 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_7 = f(t_n+h, y_{p,n+1})`
         cur_time=cur_time + dt, &
         cur_state=new_state)
      this%next_dot_state_assigned = .True.
      this%next_dot_state = dot_state_k7

      extra_state = cur_state &                                      ! `y_{p,n+1} =  y_n + \frac{5179}{57600}h k_1 + \frac{7571}{16695}h k_3 + \frac{393}{640}h k_4 - \frac{92097}{339200}h k_5 + \frac{187}{2100}h k_6 + \frac{1}{40}h k_7`
                  + dt*(5179.0_dp/57600.0_dp*dot_state_k1 + 7571.0_dp/16695.0_dp*dot_state_k3 &
                  + 393.0_dp/640.0_dp*dot_state_k4 - 92097.0_dp/339200.0_dp*dot_state_k5 &
                  + 187.0_dp/2100.0_dp*dot_state_k6 + 1.0_dp/40.0_dp*dot_state_k7)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module RK45_Dormand_Prince_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module RK45_Fehlberg_Class
   use General_Settings, only: dp, setting_max_internal_states
   use Solver_Baseclass, only: Solver
   use Constitutive_Model_Baseclass, only: Constitutive_Model
   implicit none

   ! --------------------------------------------------------------- !
   type, extends(Solver) :: RK45_Fehlberg
   ! --------------------------------------------------------------- !
      private

      contains

      procedure :: Integration_Algorithm
   end type


   contains


   ! --------------------------------------------------------------- !
   function Integration_Algorithm(this, integ_obj, cur_state, dot_state_k1, cur_time, dt) result(new_state)
   ! --------------------------------------------------------------- !
      class(RK45_Fehlberg), intent(inout) :: this
      class(Constitutive_Model), intent(inout) :: integ_obj
      real(dp), dimension(setting_max_internal_states), intent(in) :: cur_state, dot_state_k1
      real(dp), intent(in) :: cur_time, dt
      real(dp), dimension(setting_max_internal_states) :: new_state
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: dot_state_k2, dot_state_k3, dot_state_k4, &
                                                          dot_state_k5, dot_state_k6, extra_state

      dot_state_k2 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_2 = f(t_n+0.25h, y_n+0.25h k_1)`
         cur_time=cur_time + 0.25_dp*dt, &
         cur_state=cur_state + dt*0.25_dp*dot_state_k1)
      dot_state_k3 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_3 = f\left(t_n+\frac{3}{8}h, y_n+\frac{3}{32}h k_1 + \frac{9}{32}h k_2\right)`
         cur_time=cur_time + 0.375_dp*dt, &
         cur_state=cur_state + dt*(3.0_dp/32.0_dp*dot_state_k1 + 9.0_dp/32.0_dp*dot_state_k2))
      dot_state_k4 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_4 = f\left(t_n+\frac{12}{13}h, y_n+\frac{1932}{2197}h k_1 - \frac{7200}{2197}h k_2 + \frac{7296}{2197}h k_3\right)`
         cur_time=cur_time + 12.0_dp/13.0_dp*dt, &
         cur_state=cur_state + dt*(1932.0_dp/2197.0_dp*dot_state_k1 &
            - 7200.0_dp/2197.0_dp*dot_state_k2 + 7296.0_dp/2197.0_dp*dot_state_k3))
      dot_state_k5 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_5 = f\left(t_n+h, y_n+\frac{439}{216}h k_1 - 8h k_2 + \frac{3680}{513}h k_3 - \frac{845}{4104}h k_4\right)`
         cur_time=cur_time + dt, &
         cur_state=cur_state + dt*(439.0_dp/216.0_dp*dot_state_k1 - 8.0_dp*dot_state_k2 &
            + 3680.0_dp/513.0_dp*dot_state_k3 - 845.0_dp/4104.0_dp*dot_state_k4))
      dot_state_k6 = integ_obj%Get_Dot_State(ref_dt=dt, &            ! `k_6 = f\left(t_n+\frac{1}{2}h, y_n-\frac{8}{27}h k_1 + 2h k_2 - \frac{3544}{2565}h k_3 + \frac{1859}{4104}h k_4 - \frac{11}{40}h k_5\right)`
         cur_time=cur_time + 0.5_dp*dt, &
         cur_state=cur_state + dt*(-8.0_dp/27.0_dp*dot_state_k1 + 2.0_dp*dot_state_k2 &
            - 3544.0_dp/2565.0_dp*dot_state_k3 + 1859.0_dp/4104.0_dp*dot_state_k4 &
            - 11.0_dp/40.0_dp*dot_state_k5))

      new_state = cur_state + dt*(25.0_dp/216.0_dp*dot_state_k1 &    ! `y_{n+1} =  y_n + \frac{25}{216}h k_1 + \frac{1408}{2565}h k_3 + \frac{2197}{4104}h k_4 - \frac{1}{5}h k_5`
                + 1408.0_dp/2565.0_dp*dot_state_k3 + 2197.0_dp/4104.0_dp*dot_state_k4 &
                - 0.2_dp*dot_state_k5)

      extra_state = cur_state + dt*(16.0_dp/135.0_dp*dot_state_k1 &  ! `y_{p,n+1} =  y_n + \frac{16}{135}h k_1 + \frac{6656}{12825}h k_3 + \frac{28561}{56430}h k_4 - \frac{9}{50}h k_5 + \frac{2}{55}h k_6`
                + 6656.0_dp/12825.0_dp*dot_state_k3 + 28561.0_dp/56430.0_dp*dot_state_k4 &
                - 9.0_dp/50.0_dp*dot_state_k5 + 2.0_dp/55.0_dp*dot_state_k6)
      call this%Set_Error_Estimate(estimation=abs(new_state - extra_state))
   end function Integration_Algorithm
end module RK45_Fehlberg_Class


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Solver_Class
   use General_Settings, only: dp, Write_Error_And_Exit
   use Solver_Baseclass
   implicit none

   private
   public :: Select_Solver, Solver


   contains


   ! --------------------------------------------------------------- !
   subroutine Select_Solver(new_solver, name)
   ! --------------------------------------------------------------- !
      use Euler_Explicit_Class
      use Euler_Richardson_Class
      use RK23_Simpson_Class
      use RK23_Bogacki_Shampine_Class
      use RK45_Cash_Karp_Class
      use RK45_Dormand_Prince_Class
      use RK45_Fehlberg_Class
      !
      class(Solver), allocatable, intent(out) :: new_solver
      character(len=*), intent(in) :: name
      ! ------------------------------------------------------------ !
      associate(identifier => name(1:7))
         if (identifier == 'Eul-Exp') then
            allocate(Euler_Explicit::new_solver)
            call new_solver%Initialize(name='Euler-Explicit', stepsize_fixed=.True.)
         !
         else if (identifier == 'Richard') then
            allocate(Euler_Richardson::new_solver)
            ! Similar to Fellin et al. (2009): Adaptive integration of constitutive rate equations, p. 3
            call new_solver%Initialize(name='Euler-Richardson', exp_stepgrow=-0.5_dp, exp_stepshrink=-0.5_dp, &
               max_stepgrow=2.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-S 23') then
            allocate(RK23_Simpson::new_solver)
            ! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Simpson)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=4.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-BS23') then
            allocate(RK23_Bogacki_Shampine::new_solver)
            ! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Bogacki-Shampine)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=3.0_dp, min_stepshrink=0.2_dp)
         !
         else if (identifier == 'RK-CK45') then
            allocate(RK45_Cash_Karp::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Cash-Karp)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else if (identifier == 'RK-DP45') then
            allocate(RK45_Dormand_Prince::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Dormand-Prince)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else if (identifier == 'RK-F 45') then
            allocate(RK45_Fehlberg::new_solver)
            ! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Fehlberg)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)
         !
         else
            call Write_Error_And_Exit('Select_Solver: Identifier >' // identifier // '< unknown')
         end if
      end associate
   end subroutine Select_Solver
end module Solver_Class


! ==================================================================================================================== !
module Xmat_Class
   use General_Settings, only: dp, setting_len_id, setting_num_statevariables, setting_max_internal_states
   use Solver_Class, only: Solver, Select_Solver
   use Constitutive_Model_Class, only: Constitutive_Model, Select_Constitutive_Model

   implicit none

   private
   public :: Xmat, Xmat_Initialize

   ! --------------------------------------------------------------- !
   type :: Xmat
   ! --------------------------------------------------------------- !
      private
      class(Solver), allocatable :: xmat_solver
      class(Constitutive_Model), allocatable :: xmat_constitutive_model

      real(dp), dimension(setting_num_statevariables) :: statevariables
      real(dp) :: assigned_time, assigned_dt
      real(dp), dimension(6) ::  assigned_dot_strain, saved_stress
      integer :: nstates

      contains

      procedure :: Import_State
      procedure :: Calculate_Results
      procedure :: Dot_Results
   end type


   contains


   ! --------------------------------------------------------------- !
   subroutine Xmat_Initialize(xmat_obj, solver_name, constitutive_model_name, material_parameters, &
      calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      class(Xmat), allocatable, intent(out) :: xmat_obj
      character(len=*), intent(in) :: solver_name
      character(len=setting_len_id), intent(in) :: constitutive_model_name
      real(dp), dimension(:), intent(in) :: material_parameters
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      if (.not. allocated(xmat_obj)) then
         allocate(Xmat::xmat_obj)
      end if

      call Select_Solver(new_solver=xmat_obj%xmat_solver, name=solver_name)
      call Select_Constitutive_Model(new_constitutive_model=xmat_obj%xmat_constitutive_model, &
         identifier=constitutive_model_name, params=material_parameters, calculate_jacobian=calculate_jacobian, &
         firstcall=firstcall)
      call xmat_obj%xmat_solver%Set_State_Mask(statemask=xmat_obj%xmat_constitutive_model%Get_Direct_Variables_Mask())
   end subroutine Xmat_Initialize


   ! --------------------------------------------------------------- !
   subroutine Import_State(this, nstates, statevariables, stress, dot_strain, totaltime, timeincrement)
   ! --------------------------------------------------------------- !
      use General_Settings, only: setting_num_statevariables, Write_Error_And_Exit
      use Debug, only: Formatval
      !
      class(Xmat), intent(inout) :: this
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: statevariables
      real(dp), dimension(6), intent(in) :: stress, dot_strain
      real(dp), intent(in) :: totaltime, timeincrement
      ! ------------------------------------------------------------ !
      if (nstates > setting_num_statevariables) then
         call Write_Error_And_Exit('Import_State: ' // Formatval('num_states', nstates) // ' requested but only ' // &
            Formatval('setting_num_statevariables', setting_num_statevariables) // &
            ' provided - increase setting_num_statevariables?')
      end if
      this%statevariables = 0.0_dp
      this%statevariables(1:nstates) = statevariables
      this%nstates = nstates
      this%assigned_dot_strain = dot_strain
      this%assigned_time = totaltime
      this%assigned_dt = timeincrement

      call this%xmat_constitutive_model%Set_Values(overall_dt=this%assigned_dt, &
         dot_strain=this%assigned_dot_strain)

      this%saved_stress = stress
   end subroutine Import_State


   ! --------------------------------------------------------------- !
   subroutine Calculate_Results(this, exportstress, exportstates, exportjacobian, exportdt)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Is_Nan, Write_Warning, setting_epsilon, &
                                  setting_restrict_initial_substep, setting_initial_substep_scale
      use Debug, only: Formatval
      use Math_Operations, only: Nonzero_Division, Norm, Pack_States, Unpack_States
      !
      class(Xmat), intent(inout) :: this
      real(dp), dimension(6), intent(out) :: exportstress
      real(dp), dimension(this%nstates), intent(out) :: exportstates
      real(dp), dimension(6, 6), intent(out) :: exportjacobian
      real(dp), intent(out) :: exportdt
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: cur_state, new_state
      real(dp), dimension(setting_num_statevariables) :: statevariables, new_statevariables
      real(dp), dimension(6, 6) :: jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables
      integer, dimension(3) :: solver_statistics
      real(dp) :: norm_dot_strain
      logical :: normal_integration, successful_normal_integration

      jac_stress = 0.0_dp
      jac_statevariables = 0.0_dp
      solver_statistics = 0

      norm_dot_strain = Norm(this%assigned_dot_strain)
      normal_integration = ((norm_dot_strain > setting_epsilon) .and. (this%assigned_dt > setting_epsilon))
      successful_normal_integration = .False.

      if (setting_restrict_initial_substep) then
         !exportdt = min(Nonzero_Division(val=setting_initial_substep_scale, fac=norm_dot_strain), &
         !   this%assigned_dt)
         exportdt = this%assigned_dt / max(1, &
            idnint(norm_dot_strain*this%assigned_dt/setting_initial_substep_scale))
      else
         exportdt = this%assigned_dt
      end if

      ! --- Packing all cur_states in a long vector
      statevariables = this%statevariables
      cur_state = Pack_States(stress=this%saved_stress, jac_stress=jac_stress, statevariables=statevariables, &
         jac_statevariables=jac_statevariables)

      if (normal_integration) then
         ! --- Do the actual integration of the constitutive routine
         call this%xmat_solver%Integrate(integ_obj=this%xmat_constitutive_model, cur_state=cur_state, &
            cur_time=this%assigned_time, dt=this%assigned_dt, new_state=new_state, mod_dt=exportdt, &
            successful_integration=successful_normal_integration)
         solver_statistics = this%xmat_solver%Get_Statistics()
      end if

      if (successful_normal_integration) then
         ! --- Unpacking of all states out of a long vector
         call Unpack_States(input_states=new_state, stress=exportstress, jac_stress=exportjacobian, &
            statevariables=new_statevariables, jac_statevariables=jac_statevariables)
      else
         ! If either the norm of the dot strain or the timeincrement is (almost) zero or the integration fails,
         ! make one call to Get_Dot_State() of the constitutive model to get the jacobian
         new_state = this%xmat_constitutive_model%Get_Dot_State(ref_dt=this%assigned_dt, &
            cur_time=this%assigned_time, cur_state=cur_state)
         ! If the time increment is greater than zero, perform a single explicit Euler step
         if (this%assigned_dt > setting_epsilon) then
            new_state = cur_state + new_state*this%assigned_dt
         end if

         ! --- Unpacking of all states out of a long vector
         call Unpack_States(input_states=new_state, stress=exportstress, jac_stress=exportjacobian, &
            statevariables=new_statevariables, jac_statevariables=jac_statevariables)
         exportstress = this%saved_stress
         new_statevariables = this%statevariables

         ! No recommendation to change the time increment (in case integration failed and the value was changed)
         exportdt = this%assigned_dt
      end if

      exportstates = new_statevariables(1:this%nstates)
      ! NOTE: Uncomment the following to add the solver statistics to exportstates. Make sure that there are enough
      ! unused state variables because they are currently overwritten without any checks
      !exportstates(this%nstates-size(solver_statistics)+1:this%nstates) = real(solver_statistics, dp)

      exportjacobian = Nonzero_Division(val=exportjacobian, fac=this%assigned_dt)

      if (any(Is_Nan(exportstress))) then
         call Write_Warning('Calculate_Results: NaN entries found in stresses. Program will probably fail shortly' &
            // char(10) // Formatval('T', exportstress))
      end if

      if (any(Is_Nan(exportjacobian))) then
         call Write_Warning('Calculate_Results: NaN entries found in jacobian. ' &
            // 'Convergence with implicit algorithm unlikely' // char(10) // Formatval('jacobian', exportjacobian))
      end if

      if (any(Is_Nan(exportstates))) then
         call Write_Warning('Calculate_Results: NaN entries found in states. Program will probably fail shortly' &
            // char(10) // Formatval('states', exportstates))
      end if
   end subroutine Calculate_Results


   ! --------------------------------------------------------------- !
   subroutine Dot_Results(this, dotstress, dotstates, jacobian)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Is_Nan, Write_Warning
      use Debug, only: Formatval
      use Math_Operations, only: Pack_States, Unpack_States
      !
      class(Xmat), intent(inout) :: this
      real(dp), dimension(6), intent(out) :: dotstress
      real(dp), dimension(this%nstates), intent(out) :: dotstates
      real(dp), dimension(6, 6), intent(out) :: jacobian
      ! ------------------------------------------------------------ !
      real(dp), dimension(setting_max_internal_states) :: cur_state, dot_state
      real(dp), dimension(setting_num_statevariables) :: statevariables, dot_statevariables, direct_variables
      logical, dimension(setting_num_statevariables) :: direct_variables_mask
      real(dp), dimension(6, 6) :: jac_stress
      real(dp), dimension(setting_num_statevariables, 6) :: jac_statevariables

      jac_stress = 0.0_dp
      jac_statevariables = 0.0_dp

      ! --- Packing all cur_states in a long vector
      statevariables = this%statevariables
      cur_state = Pack_States(stress=this%saved_stress, jac_stress=jac_stress, statevariables=statevariables, &
         jac_statevariables=jac_statevariables)

      ! --- Calculate the response of the constitutive model
      dot_state = this%xmat_constitutive_model%Get_Dot_State(ref_dt=this%assigned_dt, cur_time=this%assigned_time, &
         cur_state=cur_state)
      call this%xmat_constitutive_model%Get_Direct_Variables(direct_variables=direct_variables)

      ! --- Unpacking of all states out of a long vector
      call Unpack_States(input_states=dot_state, stress=dotstress, jac_stress=jacobian, &
         statevariables=dot_statevariables, jac_statevariables=jac_statevariables)
      direct_variables_mask = this%xmat_constitutive_model%Get_Direct_Variables_Mask()
      where (direct_variables_mask)
         dot_statevariables =  direct_variables
      end where
      dotstates = dot_statevariables(1:this%nstates)

      if (any(Is_Nan(dotstress))) then
         call Write_Warning('Dot_Results: NaN entries found in stresses.' // char(10) // Formatval('T', dotstress))
      end if

      if (any(Is_Nan(jacobian))) then
         call Write_Warning('Dot_Results: NaN entries found in jacobian. ' &
            // 'Convergence with implicit algorithm unlikely' // char(10) // Formatval('jacobian', jacobian))
      end if

      if (any(Is_Nan(dotstates))) then
         call Write_Warning('Dot_Results: NaN entries found in states.' // char(10) // Formatval('states', dotstates))
      end if
   end subroutine Dot_Results
end module Xmat_Class


! ==================================================================================================================== !
module Custom_Utilities
   use General_Settings, only: dp, setting_len_id, setting_id_hypo_vw96, setting_id_viscohypo_ni03, &
                               setting_id_barodesy_ko21
   implicit none


   contains


   ! --------------------------------------------------------------- !
   pure function Is_First_Call(step, iteration)
   ! --------------------------------------------------------------- !
      integer, intent(in) :: step, iteration
      logical :: Is_First_Call
      ! ------------------------------------------------------------ !
      if ((step == 1) .and. (iteration == 1)) then
         Is_First_Call = .True.
      else
         Is_First_Call = .False.
      end if
   end function Is_First_Call


   ! --------------------------------------------------------------- !
   function Preprocess_State_Matrices(identifier, intername, nstates, states, rotation) result(mod_states)
   ! --------------------------------------------------------------- !
      use General_Settings, only: Write_Error_And_Exit
      use Math_Operations, only: Vec9_To_Mat, Matrix_Rotation, Mat_To_Vec9
      !
      character(len=setting_len_id), intent(in) :: identifier
      character(len=13), intent(in) :: intername
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: states
      real(dp), dimension(3, 3), intent(in) :: rotation
      real(dp), dimension(nstates) :: mod_states
      ! ------------------------------------------------------------ !
      real(dp), dimension(6) :: save_mat
      real(dp), dimension(9) :: save_vec

      mod_states = states

      if ((identifier == setting_id_hypo_vw96) .or. (identifier == setting_id_viscohypo_ni03)) then
         ! Adjust intergranular strain
         if (nstates < 10) then
            call Write_Error_And_Exit('Preprocess_State_Matrices: Expecting at least ' // &
               '10 (internal) states for this constitutive model')
         end if

         save_vec = [mod_states(2:4), 0.5_dp*mod_states(5:10)]

         ! Rotate previous intergranular strain to new coordinate orientation
         save_mat = Vec9_To_Mat(vec9=save_vec)
         save_mat = Matrix_Rotation(save_mat, rotation)
         save_vec = Mat_To_Vec9(mat=save_mat)

         if (intername == 'umat_abq_std ') then
            ! Switch order corresponding to Abaqus/Standard after transformation based on internal representation
            ! (1, 1), (2, 2), (3, 3), (1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1)
            save_vec(6:9) = [save_vec(8), save_vec(9), save_vec(6), save_vec(7)]
         end if

         mod_states(2:10) = save_vec
      else if (identifier == setting_id_barodesy_ko21) then
         ! Adjust last saved strain
         if (nstates < 12) then
            call Write_Error_And_Exit('Preprocess_State_Matrices: Expecting at least ' // &
               '12 (internal) states for this constitutive model')
         end if

         ! Rotate last saved strain to new coordinate orientation
         save_mat = Vec9_To_Mat(vec9=mod_states(4:12))
         save_mat = Matrix_Rotation(save_mat, rotation)
         mod_states(4:12) = Mat_To_Vec9(mat=save_mat)
      end if
   end function Preprocess_State_Matrices


   ! --------------------------------------------------------------- !
   pure function Postprocess_State_Matrices(identifier, intername, nstates, states) result(mod_states)
   ! --------------------------------------------------------------- !
      character(len=setting_len_id), intent(in) :: identifier
      character(len=13), intent(in) :: intername
      integer, intent(in) :: nstates
      real(dp), dimension(nstates), intent(in) :: states
      real(dp), dimension(nstates) :: mod_states
      ! ------------------------------------------------------------ !
      real(dp), dimension(9) :: ig_vec

      mod_states = states

      if ((identifier == setting_id_hypo_vw96) .or. (identifier == setting_id_viscohypo_ni03)) then
         ! Adjust intergranular strain
         ig_vec = [mod_states(2:4), 2.0_dp*mod_states(5:10)]

         if (intername == 'umat_abq_std ') then
            ! Switch order corresponding to Abaqus/Standard based on internal representation
            ! (1, 1), (2, 2), (3, 3), (1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1)
            ig_vec(6:9) = [ig_vec(8), ig_vec(9), ig_vec(6), ig_vec(7)]
         end if
         mod_states(2:10) = ig_vec
      end if
   end function Postprocess_State_Matrices
end module Custom_Utilities


#ifdef PLAXIS_DLL
! ------------------------------------------------------------------ ! ----------------------------------------------- !
module PlaxisInformationPool
   use General_Settings, only: dp, setting_num_statevariables, setting_max_mat_params, setting_len_id, &
                               setting_id_elasticity, setting_id_hypo_wu92, setting_id_hypo_vw96, &
                               setting_id_viscohypo_ni03, setting_id_barodesy_ko15, setting_id_barodesy_sc18, &
                               setting_id_barodesy_ko21
   implicit none

   ! ATTENTION: The order of the constitutive models has to be the same for all variables. Currently, if any changes
   !            are made to the constitutive models (or which constitutive models are included) adjustments may be
   !            necessary here
   character(len=setting_len_id), parameter, dimension(7) :: ModelList = [ &
      setting_id_elasticity, setting_id_hypo_wu92, setting_id_hypo_vw96, setting_id_viscohypo_ni03, &
      setting_id_barodesy_ko15, setting_id_barodesy_sc18, setting_id_barodesy_ko21]
   integer, parameter, dimension(7) :: numParameters =     [2, 4, 16, 15, 7, 9, 10]
   integer, parameter, dimension(7) :: numStateVariables = [0, 1, 11, 11, 1, 1, 12]

   character(len=16), parameter, dimension(16, 7) :: ModelParameterNamesUnits = reshape([ &
      ! Elasticity
      'E        :F/L^2#', '@nu#     :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-Wu92
      'C_1#     :-     ', 'C_2#     :-     ', 'C_3#     :-     ', 'C_4#     :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-VW96
      '@j#_c#   :rad   ', '@nu#     :-     ', 'h_s#     :F/L^2#', 'n_H      :-     ', &
      'e_d0#    :-     ', 'e_c0#    :-     ', 'e_i0#    :-     ', '@a#_H#   :-     ', &
      '@b#_H#   :-     ', 'm_T#     :-     ', 'm_R#     :-     ', 'R_max#   :-     ', &
      '@b#_r#   :-     ', '@c#      :-     ', '-        :-     ', 'e_0#     :-     ', &
      !
      ! ViHy-Ni03
      'e_100#   :-     ', '@nu#     :-     ', '@l#      :-     ', '@k#      :-     ', &
      '@b#_b#   :-     ', 'I_v#     :-     ', 'D_r#     :-     ', '@j#_c#   :rad   ', &
      '-        :-     ', 'm_T#     :-     ', 'm_R#     :-     ', 'R_max#   :-     ', &
      '@b#_r#   :-     ', '@c#      :-     ', 'OCR      :-     ', '         :-     ', &
      !
      ! Baro-Ko15
      '@j#_c#   :rad   ', 'c_2#     :-     ', 'c_3#     :-     ', 'c_4#     :-     ', &
      'c_5#     :-     ', 'e_c0#    :-     ', 'e_min#   :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Sc18
      '@j#_c#   :rad   ', 'k_r#     :-     ', '@x#      :-     ', '@k#      :-     ', &
      'e_c0#    :-     ', 'c_e#     :-     ', 'c_5#     :-     ', 'c_6#     :-     ', &
      'c_7#     :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Ko21
      '@j#_c#   :rad   ', 'c_2#     :-     ', 'c_3#     :-     ', 'c_4#     :-     ', &
      'c_5#     :-     ', '@k_1#    :-     ', '@k_2#    :-     ', 'k_1#     :-     ', &
      'k_2#     :-     ', 'e_c0#    :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ' &
      ], [16, 7])

   character(len=16), parameter, dimension(12, 7) :: ModelStatevarNamesUnits = reshape([ &
      ! Elasticity
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-Wu92
      'voidratio:-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-VW96
      'voidratio:-     ', 'igran11  :-     ', 'igran22  :-     ', 'igran33  :-     ', &
      'igran12  :-     ', 'igran21  :-     ', 'igran23  :-     ', 'igran32  :-     ', &
      'igran13  :-     ', 'igran31  :-     ', 'young_rep:F/L^2#', '         :-     ', &
      !
      ! ViHy-Ni03
      'voidratio:-     ', 'igran11  :-     ', 'igran22  :-     ', 'igran33  :-     ', &
      'igran12  :-     ', 'igran21  :-     ', 'igran23  :-     ', 'igran32  :-     ', &
      'igran13  :-     ', 'igran31  :-     ', 'young_rep:F/L^2#', '         :-     ', &
      !
      ! Baro-Ko15
      'voidratio:-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Sc18
      'voidratio:-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Ko21
      'voidratio:-     ', 'e_c1#    :-     ', 'e_c2#    :-     ', 'D_old#11 :-     ', &
      'D_old#22 :-     ', 'D_old#33 :-     ', 'D_old#12 :-     ', 'D_old#21 :-     ', &
      'D_old#23 :-     ', 'D_old#32 :-     ', 'D_old#13 :-     ', 'D_old#31 :-     ' &
      ], [12, 7])
end module PlaxisInformationPool


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
   real(dp), dimension(6) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(6, 6) :: jacobian
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
      material_parameters=materialproperties(1:nparams), calculate_jacobian=.True., firstcall=firstcall)
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


! ------------------------------------------------------------------ !
subroutine GetParamName(iModel, iParameter, ParameterName)           ! PLAXIS: Returns name of selected parameter of constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelParameterNamesUnits, numParameters
   implicit none
   !
   interface
      subroutine GetStateVarNameInternal(iModel, iStatevar, StatevarName)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarName
      end subroutine GetStateVarNameInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iParameter
   character(len=*), intent(out) :: ParameterName
   ! --------------------------------------------------------------- !
   integer :: lenName, num_params
   character(len=16) :: tempParameterName

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetParamName         ! Export function name in dll
#endif

   tempParameterName = '                '
   if ((iModel >= 1) .and. (iModel <= size(ModelParameterNamesUnits, 2))) then
      ! List not only parameters but also state variables to be able to initialize the latter by additional values in the former
      num_params = numParameters(iModel)
      if (iParameter > num_params) then
         call GetStateVarNameInternal(iModel=iModel, iStatevar=iParameter-num_params, StatevarName=tempParameterName)
         tempParameterName = tempParameterName(2:10)
         if (tempParameterName /= '') then
            tempParameterName = 'state: ' // tempParameterName
         end if
      else
         if ((iParameter >= 1) .and. (iParameter <= numParameters(iModel)) &
            .and. (iParameter <= size(ModelParameterNamesUnits, 1))) then

            ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
            tempParameterName = ModelParameterNamesUnits(iParameter, iModel)
            tempParameterName(10:16) = '       '
         end if
      end if
   end if

   lenName = len_trim(tempParameterName)
   ! The first char represents the binary length of the string
   ParameterName = char(lenName) // tempParameterName(1:lenName)
end subroutine GetParamName


! ------------------------------------------------------------------ !
subroutine GetParamUnit(iModel, iParameter, ParameterUnit)           ! PLAXIS: Returns unit of parameter of chosen constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelParameterNamesUnits, numParameters
   implicit none
   !
   interface
      subroutine GetStateVarUnitInternal(iModel, iStatevar, StatevarUnit)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarUnit
      end subroutine GetStateVarUnitInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iParameter
   character(len=*) :: ParameterUnit
   ! --------------------------------------------------------------- !
   integer :: lenUnit, num_params
   character(len=16) :: tempParameterUnit

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetParamUnit         ! Export function name in dll
#endif

   tempParameterUnit = '-               '
   if ((iModel >= 1) .and. (iModel <= size(ModelParameterNamesUnits, 2))) then
      ! List not only parameters but also state variables to be able to initialize the latter by additional values in the former
      num_params = numParameters(iModel)
      if (iParameter > num_params) then
         call GetStateVarUnitInternal(iModel=iModel, iStatevar=iParameter-num_params, StatevarUnit=tempParameterUnit)
         tempParameterUnit = tempParameterUnit(2:16) // ' '
      else
         if ((iParameter >= 1) .and. (iParameter <= numParameters(iModel)) &
            .and. (iParameter <= size(ModelParameterNamesUnits, 1))) then

            ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
            tempParameterUnit = ModelParameterNamesUnits(iParameter, iModel)
            tempParameterUnit = tempParameterUnit(11:16) // '          '
         end if
      end if
   end if
   lenUnit = len_trim(tempParameterUnit)
   ! The first char represents the binary length of the string
   ParameterUnit = char(lenUnit) // tempParameterUnit(1:lenUnit)
end subroutine GetParamUnit


! ------------------------------------------------------------------ !
subroutine GetStateVarCount(iModel, nStatevars)                      ! PLAXIS: Returns number of state variables of constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(out) :: nStatevars
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarCount     ! Export function name in dll
#endif

   nStatevars = 0
   if ((iModel >= 1) .and. (iModel <= size(numStateVariables))) then
      nStatevars = numStateVariables(iModel)
   end if
end subroutine GetStateVarCount


! ------------------------------------------------------------------ !
subroutine GetStateVarName(iModel, iStatevar, StatevarName)          ! PLAXIS: Returns name of state variable of chosen constitutive model
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   interface
      subroutine GetStateVarNameInternal(iModel, iStatevar, StatevarName)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarName
      end subroutine GetStateVarNameInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarName
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarName      ! Export function name in dll
#endif

   call GetStateVarNameInternal(iModel=iModel, iStatevar=iStatevar, StatevarName=StatevarName)
end subroutine GetStateVarName


! ------------------------------------------------------------------ !
subroutine GetStateVarNameInternal(iModel, iStatevar, StatevarName)
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarName
   ! --------------------------------------------------------------- !
   integer :: lenName
   character(len=16) :: tempStatevarName

   tempStatevarName = '                '
   if ((iModel >= 1) .and. (iModel <= size(ModelStatevarNamesUnits, 2))) then
      if ((iStatevar >= 1) .and. (iStatevar <= numStateVariables(iModel)) &
         .and. (iStatevar <= size(ModelStatevarNamesUnits, 1))) then

         ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
         tempStatevarName = ModelStatevarNamesUnits(iStatevar, iModel)
         tempStatevarName(10:16) = '       '
      end if
   end if
   lenName = len_trim(tempStatevarName)
   ! The first char represents the binary length of the string
   StatevarName = char(lenName) // tempStatevarName(1:lenName)
end subroutine GetStateVarNameInternal


! ------------------------------------------------------------------ !
subroutine GetStateVarUnit(iModel, iStatevar, StatevarUnit)          ! PLAXIS: Returns unit of state variable of constitutive model
! ------------------------------------------------------------------ !
   implicit none
   !
   interface
      subroutine GetStateVarUnitInternal(iModel, iStatevar, StatevarUnit)
         integer, intent(in) :: iModel
         integer, intent(in) :: iStatevar
         character(len=*), intent(out) :: StatevarUnit
      end subroutine GetStateVarUnitInternal
   end interface
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarUnit
   ! --------------------------------------------------------------- !
#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall,Reference :: GetStateVarUnit      ! Export function name in dll
#endif

   call GetStateVarUnitInternal(iModel=iModel, iStatevar=iStatevar, StatevarUnit=StatevarUnit)
end subroutine GetStateVarUnit


! ------------------------------------------------------------------ !
pure subroutine GetStateVarUnitInternal(iModel, iStatevar, StatevarUnit)
! ------------------------------------------------------------------ !
   use PlaxisInformationPool, only: ModelStatevarNamesUnits, numStateVariables
   implicit none
   !
   integer, intent(in) :: iModel
   integer, intent(in) :: iStatevar
   character(len=*), intent(out) :: StatevarUnit
   ! --------------------------------------------------------------- !
   integer :: lenUnit
   character(len=16) :: tempStatevarUnit

   tempStatevarUnit = '-               '
   if ((iModel >= 1) .and. (iModel <= size(ModelStatevarNamesUnits, 2))) then
      if ((iStatevar >= 1) .and. (iStatevar <= numStateVariables(iModel)) &
         .and. (iStatevar <= size(ModelStatevarNamesUnits, 1))) then

         ! The last check is just a sanity check in case the amount of elements don't match anymore after making changes
         tempStatevarUnit = ModelStatevarNamesUnits(iStatevar, iModel)
         tempStatevarUnit = tempStatevarUnit(11:16) // '          '
      end if
   end if
   lenUnit = len_trim(tempStatevarUnit)
   ! The first char represents the binary length of the string
   StatevarUnit = char(lenUnit) // tempStatevarUnit(1:lenUnit)
end subroutine GetStateVarUnitInternal
#endif


#ifdef MATLAB_CALLING
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
   real(dp), dimension(6) :: inp_stress, inp_strain, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt
   character(len=setting_len_id) :: matname

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(6, 6) :: jacobian_internal
   character(len=setting_len_id) :: identifier
   character(len=13), parameter :: intername = 'mex_function '
   integer :: idx
   logical :: firstcall

   ! The function should be called in matlab with the following seven arguments (in order):
   ! materialname, materialproperties, inp_state, inp_stress, inp_dot_strain, timeincrement, totaltime
   if (nrhs .ne. 7) then
      call Write_Error_And_Exit('mexFunction: Expecting exactly 7 input arguments from Matlab')
   end if
   ! The function returns the following three arguments to matlab (in order):
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
      material_parameters=materialproperties(1:nel_props), calculate_jacobian=.True., firstcall=firstcall)
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
#endif


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
   real(dp), dimension(6) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, totaltime, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(6, 6) :: jacobian
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
      material_parameters=materialproperties, calculate_jacobian=.True., firstcall=firstcall)
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
   real(dp), dimension(6) :: inp_stress, inp_dot_strain
   real(dp) :: timeincrement, steptime_internal, totaltime_internal, exportdt

   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(6, 6) :: jacobian
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
   real(dp), dimension(6, 6) :: jacobian_internal
   real(dp), dimension(6) :: inp_stress, inp_dot_strain, dotstress_internal
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


! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                       x m a t _ c o n s o l e                      ! Entry point for Fortran
! ------------------------------------------------------------------ !
subroutine xmat_console( &
   ! ============= variables passed in for information ============= !
      materialname, nparams, materialparameters, nstatevar, statevariables, &
      ncomponents, oldstress, oldstrain, timeincrement, totaltime, &
   ! ==================== user coding to define ==================== !
      newstress, newstate, jacobian)
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
   real(dp), dimension(ncomponents), intent(out) :: newstress
   real(dp), dimension(nstatevar), intent(out) :: newstate
   real(dp), dimension(ncomponents, ncomponents), intent(out) :: jacobian
   ! --------------------------------------------------------------- !
   class(Xmat), allocatable :: xmat_obj
   real(dp), dimension(6, 6) :: jacobian_internal
   real(dp), dimension(6) :: inp_stress, inp_dot_strain
   real(dp), dimension(ncomponents) :: tmp_dot_strain
   character(len=setting_len_id) :: identifier
   character(len=13), parameter :: intername = 'xmat_console '
   real(dp) :: exportdt
   integer :: num_dimensions, num_shear
   logical :: firstcall

   ! Determine a supported set of direct and shear components out of ncomponents
   call Estimate_Components(nel=ncomponents, num_dimensions=num_dimensions, num_shear=num_shear)
   if (.not. Check_Input_Dimensions(num_dimensions=num_dimensions, num_shear=num_shear)) then
      call Write_Error_And_Exit('xmat_console: ' // Formatval('ncomponents', ncomponents) &
         // ' can not be converted to a suitable set of direct and shear components')
   end if

   ! --- Assign to custom precision variables (consider element positioning)
   identifier = Uppercase(text=materialname(1:setting_len_id))

   inp_stress = Import_Matrix(mat=real(oldstress, dp), num_dimensions=num_dimensions, num_shear=num_shear)
   tmp_dot_strain = Nonzero_Division(val=real(oldstrain, dp), fac=timeincrement)
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
   call xmat_obj%Calculate_Results(exportstress=inp_stress, &        ! Do calculation and return values
      exportstates=newstate, exportjacobian=jacobian_internal, exportdt=exportdt)

   ! --- Reassign to calling program variables
   newstress = Export_Matrix(mat=inp_stress, num_dimensions=num_dimensions, num_shear=num_shear)
   jacobian = Export_Tensor(tens=jacobian_internal, num_dimensions=num_dimensions, num_shear=num_shear)
end subroutine xmat_console


#ifdef CINTER
! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                 D o t _ V a l u e s _ C I n t e r                  ! Dot_Values wrapper function for interoperability with C
! ------------------------------------------------------------------ !
subroutine Dot_Values_CInter( &
   ! ============= variables passed in for information ============= !
      c_materialname, c_nparams, c_materialparameters, c_nstatevar, c_statevariables, &
      c_ncomponents, c_oldstress, c_oldstrain, c_timeincrement, c_totaltime, &
   ! ==================== user coding to define ==================== !
      c_dotstress, c_dotstate, c_jacobian) bind(C, name='Dot_Values_CInter')
! ------------------------------------------------------------------ !
   use iso_c_binding, only: c_int, c_double, c_char
   use General_Settings, only: dp, C_To_Fortran_String
   implicit none

   interface
      subroutine Dot_Values(materialname, nparams, materialparameters, nstatevar, statevariables, &
         ncomponents, oldstress, oldstrain, timeincrement, totaltime, dotstress, dotstate, jacobian)
         use General_Settings, only: dp
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
      end subroutine Dot_Values
   end interface

   character(kind=c_char, len=1), dimension(80), intent(in) :: c_materialname
   integer(c_int), intent(in) :: c_nparams, c_nstatevar, c_ncomponents
   real(c_double), dimension(c_nparams), intent(in) :: c_materialparameters
   real(c_double), dimension(c_nstatevar), intent(in) ::c_statevariables
   real(c_double), dimension(c_ncomponents), intent(in) :: c_oldstress, c_oldstrain
   real(c_double), intent(in) :: c_timeincrement, c_totaltime
   real(c_double), dimension(c_ncomponents), intent(out) :: c_dotstress
   real(c_double), dimension(c_nstatevar), intent(out) :: c_dotstate
   real(c_double), dimension(c_ncomponents*c_ncomponents), intent(out) :: c_jacobian
   ! --------------------------------------------------------------- !
   character(len=80) :: materialname
   real(dp), dimension(c_nparams) :: materialparameters
   real(dp), dimension(c_nstatevar) :: statevariables, dotstate
   real(dp), dimension(c_ncomponents) :: oldstress, oldstrain, dotstress
   real(dp), dimension(c_ncomponents, c_ncomponents) :: jacobian
   real(dp) :: timeincrement, totaltime
   integer :: nparams, nstatevar, ncomponents, idx, jdx

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall :: Dot_Values_CInter              ! Export function name in dll
#endif

   materialname       = C_To_Fortran_String(length=80, c_string=c_materialname)
   nparams            = c_nparams
   materialparameters = real(c_materialparameters, dp)
   nstatevar          = c_nstatevar
   statevariables     = real(c_statevariables, dp)
   ncomponents        = c_ncomponents
   oldstress          = real(c_oldstress, dp)
   oldstrain          = real(c_oldstrain, dp)
   timeincrement      = real(c_timeincrement, dp)
   totaltime          = real(c_totaltime, dp)

   call Dot_Values(materialname=materialname, nparams=nparams, materialparameters=materialparameters, &
      nstatevar=nstatevar, statevariables=statevariables, ncomponents=ncomponents, oldstress=oldstress, &
      oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime, dotstress=dotstress, &
      dotstate=dotstate, jacobian=jacobian)

   c_dotstress = real(dotstress, c_double)
   c_dotstate = real(dotstate, c_double)
   do idx = 1, ncomponents
      do jdx = 1, ncomponents
         c_jacobian((jdx-1)*ncomponents+idx) = real(jacobian(jdx, idx), c_double)
      end do
   end do
end subroutine Dot_Values_CInter


! ------------------------------------------------------------------ ! ----------------------------------------------- !
!                x m a t _ c o n s o l e _ C I n t e r               ! xmat_console wrapper function for interoperability with C
! ------------------------------------------------------------------ !
subroutine xmat_console_CInter( &
   ! ============= variables passed in for information ============= !
      c_materialname, c_nparams, c_materialparameters, c_nstatevar, c_statevariables, &
      c_ncomponents, c_oldstress, c_oldstrain, c_timeincrement, c_totaltime, &
   ! ==================== user coding to define ==================== !
      c_newstress, c_newstate, c_jacobian) bind(C, name='xmat_console_CInter')
! ------------------------------------------------------------------ !
   use iso_c_binding, only: c_int, c_double, c_char, c_null_char
   use General_Settings, only: dp, C_To_Fortran_String
   implicit none

   interface
      subroutine xmat_console(materialname, nparams, materialparameters, nstatevar, statevariables, &
         ncomponents, oldstress, oldstrain, timeincrement, totaltime, newstress, newstate, jacobian)
         use General_Settings, only: dp
         character(len=80), intent(in) :: materialname
         integer, intent(in) :: nparams, nstatevar
         real(dp), dimension(nparams), intent(in) :: materialparameters
         real(dp), dimension(nstatevar), intent(in) :: statevariables
         integer, intent(in) :: ncomponents
         real(dp), dimension(ncomponents), intent(in) :: oldstress, oldstrain
         real(dp), intent(in) :: timeincrement, totaltime
         real(dp), dimension(ncomponents), intent(out) :: newstress
         real(dp), dimension(nstatevar), intent(out) :: newstate
         real(dp), dimension(ncomponents, ncomponents), intent(out) :: jacobian
      end subroutine xmat_console
   end interface

   character(kind=c_char, len=1), dimension(80), intent(in) :: c_materialname
   integer(c_int), intent(in) :: c_nparams, c_nstatevar, c_ncomponents
   real(c_double), dimension(c_nparams), intent(in) :: c_materialparameters
   real(c_double), dimension(c_nstatevar), intent(in) ::c_statevariables
   real(c_double), dimension(c_ncomponents), intent(in) :: c_oldstress, c_oldstrain
   real(c_double), intent(in) :: c_timeincrement, c_totaltime
   real(c_double), dimension(c_ncomponents), intent(out) :: c_newstress
   real(c_double), dimension(c_nstatevar), intent(out) :: c_newstate
   real(c_double), dimension(c_ncomponents*c_ncomponents), intent(out) :: c_jacobian
   ! --------------------------------------------------------------- !
   character(len=80) :: materialname
   real(dp), dimension(c_nparams) :: materialparameters
   real(dp), dimension(c_nstatevar) :: statevariables, newstate
   real(dp), dimension(c_ncomponents) :: oldstress, oldstrain, newstress
   real(dp), dimension(c_ncomponents, c_ncomponents) :: jacobian
   real(dp) :: timeincrement, totaltime
   integer :: nparams, nstatevar, ncomponents, idx, jdx

#ifndef NOBIB
!DEC$ ATTRIBUTES DLLExport,StdCall :: xmat_console_CInter            ! Export function name in dll
#endif

   materialname       = C_To_Fortran_String(length=80, c_string=c_materialname)
   nparams            = c_nparams
   materialparameters = real(c_materialparameters, dp)
   nstatevar          = c_nstatevar
   statevariables     = real(c_statevariables, dp)
   ncomponents        = c_ncomponents
   oldstress          = real(c_oldstress, dp)
   oldstrain          = real(c_oldstrain, dp)
   timeincrement      = real(c_timeincrement, dp)
   totaltime          = real(c_totaltime, dp)

   call xmat_console(materialname=materialname, nparams=nparams, materialparameters=materialparameters, &
      nstatevar=nstatevar, statevariables=statevariables, ncomponents=ncomponents, oldstress=oldstress, &
      oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime, newstress=newstress, &
      newstate=newstate, jacobian=jacobian)

   c_newstress = real(newstress, c_double)
   c_newstate = real(newstate, c_double)
   do idx = 1, ncomponents
      do jdx = 1, ncomponents
         c_jacobian((jdx-1)*ncomponents+idx) = real(jacobian(jdx, idx), c_double)
      end do
   end do
end subroutine xmat_console_CInter
#endif


