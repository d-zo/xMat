
xMat v0.858
===========

The [xMat](https://github.com/d-zo/xMat)
is a container structure for rate-dependent constitutive models for soil.
An object-oriented approach allows to freely combine constitutive models,
integration methods and interfaces.
In the current version it includes hypoplastic and barodetic constitutive models,
different explicit Runge-Kutta integration methods with adaptive step size control,
and interfaces to different programs and programming languages.
It was created at the Institute of Geotechnical Engineering
at the Hamburg University of Technology.



Overview
--------

The xMat is designed as a user routine and written in Fortran.
Typically, the file `xmat.f` is used directly or compiled in some way with other programs.
It was originally designed to calculate the stress response of a soil element given some specific initial conditions in simple Finite-Element-Simulations.

There is no warranty for the xMat and the correctness of the calculated results (see also `LICENSE` file).
But if you encounter a bug or flawed behaviour, please report it if it wasn't reported yet.



Using the xMat
--------------

The xMat is a library which needs a driver program calling one of the interfaces.
Then depending on the program, operating system, and compiler to use, it can be used in various ways.
Some examples, wrapper scripts, and notes for different interfaces are provided in `tests/interface-tests/input/`
(including usage of the xMat with Plaxis or Abaqus).

But before compiling and using the xMat, it is recommended to have a look at the settings in `xmat.f` and possibly adjust them.
All relevant internal parameter are defined in the module `General_Settings` so that they can be adjusted at one central place.

More documentation (in German) can be found in the PhD thesis

 - D. Zobel (2023): Zur Implementierung höherwertiger Stoffmodelle am Beispiel der Hypoplastizität
   (link to document will follow after publication)



Working with xMat in Linux
--------------------------

Linux is the main development platform for xMat.
To work properly, Unix Makefiles and Python 3.5+ have to be supported.
The `Makefile` can execute (among others) the following tasks:

 - `make prepare` to select the internal matrix and tensor representation and create necessary symlinks
 - `make xmat` to assemble `xmat.f` from source files
 - `make examples` to compile an example program with GNU Fortran compiler (`gfortran`)
   and Intel Fortran compiler (`ifort`). This is also the default when calling `make` without an argument,
   i.e. assemble `xmat.f` and compile example programs
 - `make test` or `make unit` to compile a program with unit tests for various xMat functions
 - `make doc` to create a LaTeX overview `xmat.pdf` and Doxygen html-documentation.
   The LaTeX overview can also be called with `make pdf` and needs a LaTeX distribution with working `pdflatex`,
   the Doxygen html-documentation with `make doxygen` and needs Doxygen


There are also some interface tests in `tests/interface-tests` which can be run in Linux and Windows.
A Fortran compiler is expected to be available and by default the GNU Fortran compiler `gfortran` is used.
This can be adjusted in the `Linux` part of the variable `settings` in `preparation.py`.
A part of this target assumes that GNU Octave is installed.
If this is not the case, remove the corresponding line in the `Makefile` before running `make`.
To create all interface tests and run them in Linux,
enter `make` and then `make run` in that folder.



Preparing xMat in Windows
-------------------------

Only a subset of the available options are also provided for developing the xMat in Windows.
Currently, using Windows the components of the xMat can be selected and assembled.
Compilation or creation of documentation is not supported and but limited testing is.


**Prerequisites**

 - Python 3.5+ installed
 - Set variable `python_path` in file `Preparation.bat` to path for actual Python executable
   (typically similar to `C:\Users\<username>\AppData\Local\Programs\Python\Python<version>\python.exe`)
 - Execute `Preparation.bat`


**Interface testing**

Below `tests/interface-tests` the files `01_Prepare.bat` and `preparation.py` have to be adjusted
manually for interface testing.
Specifically, the path to the Python interpreter has to be inserted in the second line of `01_Prepare.bat`
and the whole `Windows`-block of the `settings`-variable in `preparation.py` has to be adjusted as well.
The default configuration assumes the Intel Fortran compiler and a path similar to an Intel oneAPI installation.



Code testing
------------

Basic code testing is done with [FRUIT](https://sourceforge.net/projects/fortranxunit/), the Fortran Unit Testing Framework.
A copy of FRUIT's file `fruit.f90` can be found in the `src/` folder of the archive `fruit_3.4.3.zip`
has to be copied into `tests/unit-tests`.
Most math functions and basic functions of the xMat are covered by the unit test framework
(the tests are located in the same directory as FRUIT).

_Note: Some test functions only have a few checks and would benefit from additional tests._

_Note: Currently the first test in `Inverse_Tensor_Test()` fails for `REPR_MAT66` - which has to be investigated in more detail.
Possibly the wrapper subroutine `Inverse_Tensor()` in xMat is not consistent with the rest._



Contributing
------------

**Bug reports**

If you found a bug, make sure you can reproduce it with the latest version of xMat.
Please check that the expected results can actually be achieved by other means
and are not considered invalid operations due to problematic template files.
Please give detailed and reproducible instructions in your report including

 - the xMat version
 - the expected result
 - the result you received
 - the command(s) used as a _minimal working example_

Note: The bug should ideally be reproducible by the _minimal working example_ alone.
Please keep the example code as short as possible (minimal).


**Feature requests**

If you have an idea for a new feature, consider searching the
[open issues](https://github.com/d-zo/xMat/issues) and
[closed issues](https://github.com/d-zo/xMat/issues?q=is%3Aissue+is%3Aclosed) first.
Afterwards, please submit a report in the
[Issue tracker](https://github.com/d-zo/xMat/issues) explaining the feature and especially

 - why this feature would be useful (use cases)
 - what could possible drawbacks be (e.g. compatibility, dependencies, ...)



License
-------

xMat is released under the
[GPL](https://www.gnu.org/licenses/gpl-3.0.html "GNU General Public License"),
version 3 or greater (see als [LICENSE](https://github.com/d-zo/xMat/blob/master/LICENSE) file).
It is provided without any warranty.

