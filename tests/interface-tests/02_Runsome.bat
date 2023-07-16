cd output

cd C
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..

cd Fortran
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..

cd Fortran-UMAT
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..

cd Fortran-USER_MOD
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..

cd Fortran-VUMAT
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..

cd Python
cmd /c 01_compile.bat nopause
cmd /c 02_run.bat nopause
cd ..
pause
