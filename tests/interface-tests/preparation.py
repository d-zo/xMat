#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Adjust paths and commands to your system
settings = {
    'Windows': {
        'Precommand': r'@call "C:\Path\to\Intel\oneAPI\compiler\version\env\vars.bat" intel64 vs2019',
        'Python': r'C:\Path\to\Python\python.exe',
        'Compiler': 'ifort',
        'Flags': '/free /fpp',
        'Precompile': '/c',
        'Dyn. library': '/dll /recursive /MT'
    },
    'Linux': {
        'Precommand': '',
        'Python': 'python',
        'Compiler': 'gfortran',
        'Flags': '-ffree-form -cpp',
        'Precompile': '-c',
        'Dyn. library': '-shared -fPIC'
    }
}


# -------------------------------------------------------------------------------------------------
def Create_Windows_Testfiles(output_dir, prog_lang, testfile, testsettings):
    command = testsettings['Compiler'] + ' ' + testsettings['Flags']
    precommand = testsettings['Precommand']

    if (prog_lang in ['Matlab', 'Octave', 'Simulink']):
        return

    with open(output_dir + '01_Compile.bat', 'w', encoding='utf-8') as outfile:
        if (prog_lang == 'C'):
            Precompile = testsettings['Precompile']
            outfile.write('''@echo off
''' + precommand + '''
''' + command + ''' /DCINTER /DNOBIB xmat.f ''' + Precompile + '''
cl xmat.obj ''' + testfile + '''.c /link /out:''' + testfile + '''.exe
del *.mod *.obj *.exp *.lib
if /i not "%~1"=="nopause" pause
''')
        elif (prog_lang == 'Fortran-USER_MOD'):
            outfile.write('''@echo off
''' + precommand + '''
''' + command + ''' /DPLAXIS_DLL /DNOBIB xmat.f ''' + testfile + '''.f /exe:''' + testfile + '''.exe
del *.mod *.obj *.exp *.lib
if /i not "%~1"=="nopause" pause
''')
        elif (prog_lang.startswith('Fortran')):
            outfile.write('''@echo off
''' + precommand + '''
''' + command + ''' xmat.f ''' + testfile + '''.f /exe:''' + testfile + '''.exe
del *.mod *.obj *.exp *.lib
if /i not "%~1"=="nopause" pause
''')
        elif (prog_lang == 'Python'):
            dyn_library = testsettings['Dyn. library']
            outfile.write('''@echo off
''' + precommand + '''
''' + command + ''' ''' + dyn_library + ''' /DCINTER xmat.f /exe:xmatc.dll
del *.mod
del xmatc.exp xmatc.lib xmat.obj
if /i not "%~1"=="nopause" pause
''')

    with open(output_dir + '02_Run.bat', 'w', encoding='utf-8') as outfile:
        if ((prog_lang == 'C') or (prog_lang.startswith('Fortran'))):
            outfile.write('''@echo off
''' + testfile + '''.exe > output.win
if /i not "%~1"=="nopause" pause
''')
        elif (prog_lang == 'Python'):
            outfile.write('''@echo off
''' + testsettings['Python'] + ''' ''' + testfile + '''.py > output.win
if /i not "%~1"=="nopause" pause
''')


# -------------------------------------------------------------------------------------------------
def Create_Linux_Testfiles(output_dir, prog_lang, testfile, testsettings):
    command = testsettings['Compiler'] + ' ' + testsettings['Flags']

    if ((prog_lang not in ['C', 'Python']) and (not prog_lang.startswith('Fortran'))):
        return

    with open(output_dir + 'Makefile', 'w', encoding='utf-8') as outfile:
        if (prog_lang == 'C'):
            Precompile = testsettings['Precompile']
            outfile.write('''.PHONY: all run

all: ''' + testfile + ''' run
            
''' + testfile + ''': xmat.f ''' + testfile + '''.c
	''' + command + ''' -DCINTER -DNOBIB $< ''' + Precompile + ''';
	gcc $(^:.f=.o) -o $@ -lgfortran -lm;
	-rm *.mod xmat.o;

run:
	./''' + testfile + ''' > output.lnx;
''')
        elif (prog_lang == 'Fortran-USER_MOD'):
            outfile.write('''.PHONY: all run

all: ''' + testfile + ''' run

''' + testfile + ''': xmat.f ''' + testfile + '''.f
	''' + command + ''' -DPLAXIS_DLL -DNOBIB $^ -o $@;
	-rm *.mod;

run:
	./''' + testfile + ''' > output.lnx;
''')
        elif (prog_lang.startswith('Fortran')):
            outfile.write('''.PHONY: all run

all: ''' + testfile + ''' run

''' + testfile + ''': xmat.f ''' + testfile + '''.f
	''' + command + ''' $^ -o $@;
	-rm *.mod;

run:
	./''' + testfile + ''' > output.lnx;
''')
        elif (prog_lang == 'Python'):
            dyn_library = testsettings['Dyn. library']
            outfile.write('''.PHONY: all run

all: xmatc.so run

xmatc.so: xmat.f
	''' + command + ''' ''' + dyn_library + ''' -DCINTER $^ -o $@;
	-rm *.mod;

run:
	''' + testsettings['Python'] + ''' ''' + testfile + '''.py > output.lnx;
''')


# -------------------------------------------------------------------------------------------------
def Create_Additional_Files(output_dir, prog_lang, testfile, system_name, testsettings):
    if (prog_lang == 'Matlab'):
        with open(output_dir + 'compile_' + system_name + '.m', 'w', encoding='utf-8') as outfile:
            outfile.write('''% Run this script in Matlab to compile xmat.f as a MEX function
% It might be necessary to configure the ''' + testsettings['Compiler'] + ''' compiler in
% Matlab first or adjust the following compilation command

%% Compilation
% This is only needed before the first use of the MEX function
% and after xmat.f has changed
mex \'COMPFLAGS="$COMPFLAGS ''' + testsettings['Flags'] + '''"\' -largeArrayDims -DMATLAB_CALLING xmat.f
delete *.mod;
''')

    if (prog_lang == 'Octave'):
        with open(output_dir + 'compile_' + system_name + '.m', 'w', encoding='utf-8') as outfile:
            outfile.write('''% Run this script in Matlab to compile xmat.f as a MEX function
% It might be necessary to configure the ''' + testsettings['Compiler'] + ''' compiler in
% Matlab first or adjust the following compilation command

%% Compilation
% This is only needed before the first use of the MEX function
% and after xmat.f has changed
mkoctfile ''' + testsettings['Flags'] + ''' -DOCTAVE_CALLING xmat.f -c;
delete *.mod;
mkoctfile -L/usr/lib64/octave/6.4.0 -o xmat_oct xmat.o xmat_oct.cpp;
''')

    elif (prog_lang == 'Simulink'):
        with open(output_dir + 'compile_' + system_name + '.m', 'w', encoding='utf-8') as outfile:
            outfile.write('''% Run this script in Matlab to compile xmat.f as a C MEX S-function
% It might be necessary to configure the ''' + testsettings['Compiler'] + ''' (and a C compiler)
% in Matlab first or adjust the following compilation commands

%% Compilation
% This is only needed before the first use of the C MEX S-function
% and after xmat.f or xmat_sl.c have changed
mex \'COMPFLAGS="$COMPFLAGS ''' + testsettings['Flags'] + '''"\' -largeArrayDims -DMATLAB_CALLING -DCINTER xmat.f -c
mex -largeArrayDims xmat.obj xmat_sl.c
delete *.mod xmat.obj;
''')


# -------------------------------------------------------------------------------------------------
def Create_Folder_If_Not_Present(foldername):
    import os

    if (not os.path.isdir(foldername)):
        os.makedirs(foldername)


# -------------------------------------------------------------------------------------------------
def main(settings):
    import os
    import glob
    import shutil

    input_folder = 'input/'
    if (not os.path.isdir(input_folder)):
        print('Error: Input folder ' + input_folder + ' does not exist')
        return

    output_folder = 'output/'
    Create_Folder_If_Not_Present(foldername=output_folder)

    testfile = 'tester'

    for cur_folder in os.listdir(input_folder):
        subfolder = input_folder + os.sep + cur_folder
        if (not os.path.isdir(subfolder)):
            continue

        output_dir = output_folder + cur_folder + os.sep
        Create_Folder_If_Not_Present(foldername=output_dir)

        relevant_files = glob.glob(subfolder + os.sep + '*')

        shutil.copy('../../xmat.f', output_dir)
        for single_file in relevant_files:
            shutil.copy(single_file, output_dir)

        # Windows files
        testsettings = settings['Windows']
        Create_Windows_Testfiles(output_dir=output_dir, prog_lang=cur_folder,
            testfile=testfile, testsettings=testsettings)

        Create_Additional_Files(output_dir=output_dir, prog_lang=cur_folder,
            testfile=testfile, system_name='win', testsettings=testsettings)

        # Linux files
        testsettings = settings['Linux']
        Create_Linux_Testfiles(output_dir=output_dir, prog_lang=cur_folder,
            testfile=testfile, testsettings=testsettings)

        Create_Additional_Files(output_dir=output_dir, prog_lang=cur_folder,
            testfile=testfile, system_name='lnx', testsettings=testsettings)


main(settings=settings)
