@echo off
set python_path="C:\Path\to\Python\python.exe"
%python_path% tools\component_selection.py
call tools\copy_common_rep_files.bat
%python_path% tools\unify_xmat.py
%python_path% tools\adjust_xmat.py xmat.f
pause
