!                                                  xMat constitutive model container v--version--
!  External program                                         D. Zobel (--cDatum--)                                         External program
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
!    -> If the reference is a xMat function, possibly the respective preprocessor macro is not defined ('-D<Makro>')
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


! Using preprocessor directives to fit routine to calling program. Abaqus and but Matlab have to explicitly set one
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


#adddirectory module General_Settings


#adddirectory module Debug


#adddirectory module Math_Operations


#adddirectory module Constitutive_Model


#adddirectory module Solver


#adddirectory module Xmat


#adddirectory module Custom_Utilities


#adddirectory hub Interfaces
