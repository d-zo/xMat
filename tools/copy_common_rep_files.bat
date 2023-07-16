@echo off
setlocal enabledelayedexpansion
pushd "%~dp0"
pushd ..\src\Math_Operations\
echo %cd%
for %%s in (common\*) do (
   echo %%s
   set rmfile=%%~nxs
   for %%d in (mat66 mat99 tens3333) do (
      if exist "%cd%\%%d\!rmfile!" del "%cd%\%%d\!rmfile!"
rem Creating symlinks seems not to work for non-admin users
rem      mklink "%cd%\%%d\!rmfile!" "%cd%\common\!rmfile!"
      copy "%cd%\common\!rmfile!" "%cd%\%%d\"
   )
)
popd
pause