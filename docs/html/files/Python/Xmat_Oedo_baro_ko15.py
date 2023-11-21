#!/usr/bin/python3
# -*- coding: utf-8 -*-

from xmat_py import xMat


maxiter = 50000

materialname = 'Baro-Ko15_Hostun Sand 01'
voidratio = 0.7
materialparameters = [0.5899, 1.0, -2.5076, 1.0, 40.0, 0.91, 0.5]
statevariables = [voidratio, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

oedo_pressure = [-10.0, -100.0, -50.0, -200.0]
timeincrement = 0.0001
totaltime = 0.0001

stress = [oedo_pressure[0]/2, oedo_pressure[0], oedo_pressure[0]/2, 0.0, 0.0, 0.0]
strain = [timeincrement*x for x in [0.0, -0.05, 0.0, 0.0, 0.0, 0.0]]

oldstress = stress
oldstate = statevariables

idx = 1
breakall = False
numbersign = 1.0


import os

if (os.name == 'posix'):
    xmat_obj = xMat(lib_path='./xmatc.so')
elif (os.name == 'nt'):
    xmat_obj = xMat(lib_path='./xmatc.dll')
else:
    print('OS not supported')
    xmat_obj = None


from timeit import default_timer as timer

with open('Xmat_Oedo_Python.csv', 'w') as outfile:
    outfile.write('{:13.6f}   {:13.6f}\n'.format(-oldstress[1], voidratio))

    start_time = timer()
    for istep in range(1, len(oedo_pressure)):
        numbersign = (-1.0)**(istep+1)
        oldstrain = [numbersign*x for x in strain]
        while (True):
            if (numbersign*oldstress[1] < numbersign*oedo_pressure[istep]):
                break

            newstress, newstate, jacobian = xmat_obj.xmat_console_CInter(materialname=materialname,
                materialparameters=materialparameters, statevariables=oldstate, oldstress=oldstress,
                oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime)

            voidratio = newstate[0]
            oldstress = newstress
            oldstate = newstate

            outfile.write('{:13.6f}   {:13.6f}\n'.format(-oldstress[1], voidratio))

            idx += 1;
            if (idx > maxiter):
                breakall = True
                break

        if (breakall):
            break

    end_time = timer()
    print('Xmat_Oedo_Python: {:6.3f}s'.format(end_time - start_time))
