#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
component_selection.py   v0.6
2020-2023 Dominik Zobel
"""


multioptioncollection = [
    ['Which Constitutive Models to include?', 'Constitutive_Model',
        [['Hypoplasticity (Wu 1992)',           'Hypoplasticity_Wu92',      False],
        ['Hypoplasticity (Niemunis 1997)',      'Hypoplasticity_VW96',      True],
        ['Viscohypoplasticity (Niemunis 2003)', 'Viscohypoplasticity_Ni03', False],
        ['Barodesy (Kolymbas 2015)',            'Barodesy_Ko15',            False],
        ['Barodesy (Schranz 2018)',             'Barodesy_Sc18',            False],
        ['Barodesy (Kolymbas 2021)',            'Barodesy_Ko21',            False],
        ['Test mode (differential equations)',  'Test_DGL',                 False]]],
    ['Which Solvers to include?', 'Solver',
        [['Euler Explicit',                    'Euler_Explicit',        True],
        ['Euler mit Richardson',               'Euler_Richardson',      False],
        ['Runge-Kutta 2/3 (Simpson)',          'RK23_Simpson',          False],
        ['Runge-Kutta 2/3 (Bogacki-Shampine)', 'RK23_Bogacki_Shampine', False],
        ['Runge-Kutta 4/5 (Cash-Karp)',        'RK45_Cash_Karp',        False],
        ['Runge-Kutta 4/5 (Dormand-Prince)',   'RK45_Dormand_Prince',   False],
        ['Runge-Kutta 4/5 (Fehlberg)',         'RK45_Fehlberg',         False]]],
    ['Which Interfaces to include?', 'Interfaces',
        [['Plaxis',         'Plaxis',          False],
        ['Matlab',          'Matlab',          False],
        ['Abaqus/Standard', 'Abaqus/Standard', True],
        ['Abaqus/Explicit', 'Abaqus/Explicit', True],
        ['Fortran',         'Fortran',         True],
        ['C',               'C',               False]]]]


math_operations_module = 'src/Math_Operations/_module.f'
singleoptioncollection = [
    'Which Representation to use?', 'Representation',
        [['Full (3x3x3x3 tensors, 3x3 matrices)',            'tens3333', False],
        ['Direct Projection (9x9 tensors, 9x1 matrices)',    'mat99',    False],
        ['Symmetric Projection (6x6 tensors, 6x1 matrices)', 'mat66',    True]]]


constitutive_model_module = 'src/Constitutive_Model/_module.f'
constitutive_model_selection = 'src/Constitutive_Model/Select_Constitutive_Model.f'
cm_names = {
    'Hypoplasticity_Wu92':      'setting_id_hypo_wu92',
    'Hypoplasticity_VW96':      'setting_id_hypo_vw96',
    'Viscohypoplasticity_Ni03': 'setting_id_viscohypo_ni03',
    'Barodesy_Ko15':            'setting_id_barodesy_ko15',
    'Barodesy_Sc18':            'setting_id_barodesy_sc18',
    'Barodesy_Ko21':            'setting_id_barodesy_ko21',
    'Test_DGL':                 'setting_id_test_dgl'}


solver_module = 'src/Solver/_module.f'
solver_selection = 'src/Solver/Select_Solver.f'
solv_options = {
    'Euler_Explicit': ['Eul-Exp',
         '''call new_solver%Initialize(name='Euler-Explicit', stepsize_fixed=.True.)'''],
    'Euler_Richardson': ['Richard',
         '''! Similar to Fellin et al. (2009): Adaptive integration of constitutive rate equations, p. 3
            call new_solver%Initialize(name='Euler-Richardson', exp_stepgrow=-0.5_dp, exp_stepshrink=-0.5_dp, &
               max_stepgrow=2.0_dp, min_stepshrink=0.2_dp)'''],
    'RK23_Simpson': ['RK-S 23',
         '''! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Simpson)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=4.0_dp, min_stepshrink=0.2_dp)'''],
    'RK23_Bogacki_Shampine': ['RK-BS23',
         '''! Inspired by from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK2/3 (Bogacki-Shampine)', exp_stepgrow=-1.0_dp/3.0_dp, &
               exp_stepshrink=-0.5_dp, max_stepgrow=3.0_dp, min_stepshrink=0.2_dp)'''],
    'RK45_Fehlberg': ['RK-F 45',
         '''! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Fehlberg)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)'''],
    'RK45_Dormand_Prince': ['RK-DP45',
         '''! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Dormand-Prince)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)'''],
    'RK45_Cash_Karp': ['RK-CK45',
         '''! Taken from Press (1997): Numerical Recipes in Fortran, p. 712
            call new_solver%Initialize(name='RK4/5 (Cash-Karp)', exp_stepgrow=-0.2_dp, &
               exp_stepshrink=-0.25_dp, max_stepgrow=5.0_dp, min_stepshrink=0.1_dp)''']}


interface_hub = 'src/Interfaces/_hub.f'
interface_includes = {
    'Abaqus/Standard': '#addfile subroutine Abaqus/UMAT\n\n\n',
    'Abaqus/Explicit': '#addfile subroutine Abaqus/VUMAT\n\n\n',
    'Matlab':          '#ifdef MATLAB_CALLING\n#addfile subroutine mexFunction\n#endif\n\n\n',
    'Plaxis':          '#ifdef PLAXIS_DLL\n#adddirectory hub Plaxis\n#endif\n\n\n',
    'Fortran':         '#addfile subroutine dot_values\n\n\n#addfile subroutine xmat_console\n\n\n',
    'C':               '#ifdef CINTER\n#addfile subroutine dot_values_cinter\n\n\n#addfile subroutine xmat_console_cinter\n#endif\n\n\n'}


# -------------------------------------------------------------------------------------------------
def Adjust_Constitutive_Models(selected_options, avail_options):
    with open(constitutive_model_module, 'w') as outfile:
        outfile.write('''\
#adddirectory module Baseclass


#adddirectory module Elasticity


''' + '\n\n\n'.join(['#adddirectory module ' + cm for cm in selected_options]) + '''


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


#addfile subroutine Select_Constitutive_Model
end module Constitutive_Model_Class
''')

    selectorblock = ''
    for cmodel in selected_options:
        selectorblock += '      else if (identifier == ' + avail_options[cmodel] + ') then\n' \
            + '         allocate(' + cmodel + '::new_constitutive_model)\n'

    with open(constitutive_model_selection, 'w') as outfile:
        outfile.write('''\
   ! --------------------------------------------------------------- !
   subroutine Select_Constitutive_Model(new_constitutive_model, identifier, params, calculate_jacobian, firstcall)
   ! --------------------------------------------------------------- !
      use Elasticity_Class
''' + '\n'.join(['      use ' + cm + '_Class' for cm in selected_options]) + '''
      !
      class(Constitutive_Model), allocatable, intent(out) :: new_constitutive_model
      character(len=setting_len_id), intent(in) :: identifier
      real(dp), dimension(:), intent(in) :: params
      logical, intent(in) :: calculate_jacobian, firstcall
      ! ------------------------------------------------------------ !
      if (identifier == setting_id_elasticity) then
         allocate(Elasticity::new_constitutive_model)
''' + selectorblock + '''\
      else
         call Write_Error_And_Exit('Select_Constitutive_Model: Identifier >' // identifier // '< unknown')
      end if
      call new_constitutive_model%Initialize(params=params, calculate_jacobian=calculate_jacobian, &
         firstcall=firstcall)
   end subroutine Select_Constitutive_Model
''')


# -------------------------------------------------------------------------------------------------
def Adjust_Solver(selected_options, avail_options):
    with open(solver_module, 'w') as outfile:
        outfile.write('''\
#adddirectory module Baseclass


''' + '\n\n\n'.join(['#adddirectory module ' + solv for solv in selected_options]) + '''


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Solver_Class
   use General_Settings, only: dp, Write_Error_And_Exit
   use Solver_Baseclass
   implicit none

   private
   public :: Select_Solver, Solver


   contains


#addfile subroutine Select_Solver
end module Solver_Class
''')

    selectorblock = ''
    for solv in selected_options:
        selectorblock += 'if (identifier == \'' + avail_options[solv][0] + '\') then\n' \
            + '            allocate(' + solv + '::new_solver)\n            ' + avail_options[solv][1] \
            + '\n         !\n         else '

    with open(solver_selection, 'w') as outfile:
        outfile.write('''\
   ! --------------------------------------------------------------- !
   subroutine Select_Solver(new_solver, name)
   ! --------------------------------------------------------------- !
''' + '\n'.join(['      use ' + solv + '_Class' for solv in selected_options]) + '''
      !
      class(Solver), allocatable, intent(out) :: new_solver
      character(len=*), intent(in) :: name
      ! ------------------------------------------------------------ !
      associate(identifier => name(1:7))
         ''' + selectorblock[:-1] + '''
            call Write_Error_And_Exit('Select_Solver: Identifier >' // identifier // '< unknown')
         end if
      end associate
   end subroutine Select_Solver
''')


# -------------------------------------------------------------------------------------------------
def Adjust_Interfaces(selected_options, avail_options):
    if (('C' in selected_options) and ('Fortran' not in selected_options)):
        print('C interfaces depend on existing Fortran interfaces - activating Fortran interfaces')
        idx_c = selected_options.index('C')
        selected_options = selected_options[:idx_c] + ['Fortran'] + selected_options[idx_c:]

    with open(interface_hub, 'w') as outfile:
        for iblock in selected_options:
            outfile.write(avail_options[iblock])


# -------------------------------------------------------------------------------------------------
def Adjust_Representation(representation):
    mod_content = '#adddirectory module ' + representation + '\n'
    aktueller_inhalt = ''
    with open(math_operations_module, 'r', encoding='utf-8') as infile:
        for zeile in infile:
            aktueller_inhalt += zeile

    if (mod_content != aktueller_inhalt):
        with open(math_operations_module, 'w', encoding='utf-8') as outfile:
            outfile.write(mod_content)


# -------------------------------------------------------------------------------------------------
def Selection_List(description, options):
    """Creates an interactive selection list with the given description. options must be a list with
    up to nine tuples of exactly two entries. The first entry will be the shown option, the second
    entry the corresponding return value if this option is selected. All options will be numbered
    starting from 1 and can be selected with their assigned value (1-9). If 0 (Exit) is chosen,
    None ist returned.
    """
    while True:
        print(description)
        for idx, singleoption in enumerate(options):
            print('      (' + str(idx+1) + ') ' + singleoption[0])

        print('      (0) Exit')
        try:
            inputvalue = int(input('> '))
        except:
            # Input was not a valid integer
            inputvalue = len(options)+1

        if (inputvalue == 0):
            return None
        elif (inputvalue > 0) and (inputvalue <= len(options)):
            return options[inputvalue-1][1]
        else:
            print('Invalid value\n')


# -------------------------------------------------------------------------------------------------
def Multi_Selection_List(description, options):
    """Creates an interactive selection list with the given description. options must be a list with
    up to nine tuples of exactly three entries. The first entry will be the shown option, the second
    entry the corresponding return value if this option is selected and the third entry a boolean
    value, if this value is preselected. All options will be numbered starting from 1 and can be
    selected with their assigned value (1-9). If 0 (Done) is chosen, a list of the return values
    of all chosen options will be returned.
    """
    while True:
        print(description)
        selectionlist = []
        for idx, singleoption in enumerate(options):
            if (singleoption[2]):
                indentchar = '    * ('
                selectionlist += [singleoption[1]]
            else:
                indentchar = '      ('

            print(indentchar + str(idx+1) + ') ' + singleoption[0])

        print('      (0) Done')
        try:
            inputvalue = int(input('> '))
        except:
            # Input was not a valid integer
            inputvalue = len(options)+1

        if (inputvalue == 0):
            return selectionlist
        elif (inputvalue > 0) and (inputvalue <= len(options)):
            # Toggle selected option
            options[inputvalue-1][2] = not options[inputvalue-1][2]
        else:
            print('Invalid value\n')


# -------------------------------------------------------------------------------------------------
def main(multioptioncollection, singleoptioncollection, cm_names, solv_options, interface_includes):
    """Expects two list-like structures with options to choose from
    """
    alloptions = dict()

    # 1) Query which representation to choose based on options in singleoptioncollection
    representation = Selection_List(description=singleoptioncollection[0], options=singleoptioncollection[2])
    if (representation is None):
        print('Exiting. Nothing changed')
        return

    Adjust_Representation(representation=representation)

    # 2) Query if single options from multioptioncollection should be selected or everything
    startoptions = [
        'Include everything?', 'select_all',
            [['yes', True, False],
            ['no, manually select components', False, True]]]
    checkall = Selection_List(description=startoptions[0], options=startoptions[2])

    if (checkall is None):
        print('Exiting. Nothing changed')
        return

    if (checkall):
        # 2a) Include everything
        kurzfassung = [[y[1] for y in x[2]] for x in multioptioncollection]
        for question, category, options in multioptioncollection:
            alloptions.update([(category, [x[1] for x in options])])
    else:
        # 2b) Query for each option in multioptioncollection which should be included
        for option in multioptioncollection:
            selection = Multi_Selection_List(description=option[0], options=option[2])
            alloptions.update([(option[1], selection)])

    # 3) Do the actual adjustments
    Adjust_Constitutive_Models(selected_options=alloptions['Constitutive_Model'], avail_options=cm_names)
    Adjust_Solver(selected_options=alloptions['Solver'], avail_options=solv_options)
    Adjust_Interfaces(selected_options=alloptions['Interfaces'], avail_options=interface_includes)
    print('Adjustments done')


main(multioptioncollection=multioptioncollection, singleoptioncollection=singleoptioncollection,
    cm_names=cm_names, solv_options=solv_options, interface_includes=interface_includes)
