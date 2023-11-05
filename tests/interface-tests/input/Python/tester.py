#!/usr/bin/python3
# -*- coding: utf-8 -*-

from xmat_py import xMat


materialname = 'HYPO-VW96_Test'
materialparameters = [0.5777, 0.0, 4000000.0, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
                      1.1, 2.2, 0.0001, 0.1, 5.5, 0.0, 0.8]
statevariables = [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

oldstress = [-100.0, -250.0, -250.0, 50.0, 10.0, 0.0]

oldstrain = [-0.00005, 0.000015, 0.000015, 0.0, 0.0, 0.0]
timeincrement = 0.001
totaltime = 0.0

import os

if (os.name == 'posix'):
   xmat_obj = xMat(lib_path='./xmatc.so')
elif (os.name == 'nt'):
   xmat_obj = xMat(lib_path='./xmatc.dll')
else:
   print('OS not supported')
   xmat_obj = None


dotstress, dotstate, jacobian = xmat_obj.Dot_Values_CInter(materialname=materialname,
   materialparameters=materialparameters, statevariables=statevariables, oldstress=oldstress,
   oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime)

print('dotstress = ' + ''.join(['{:16.7f}'.format(x) for x in dotstress]))
print('dotstate = ' + ''.join(['{:16.7f}'.format(x) for x in dotstate]))

print('jacobian = ')
for idx in range(6):
   print('( ' + str(idx+1) + ', :) ' + ''.join(['{:16.7f}'.format(x) for x in jacobian[6*idx:6*(idx+1)]]))

print('')

newstress, newstate, jacobian = xmat_obj.xmat_console_CInter(materialname=materialname,
   materialparameters=materialparameters, statevariables=statevariables, oldstress=oldstress,
   oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime)

print('newstress = ' + ''.join(['{:16.7f}'.format(x) for x in newstress]))
print('newstate = ' + ''.join(['{:16.7f}'.format(x) for x in newstate]))

print('jacobian = ')
for idx in range(6):
   print('( ' + str(idx+1) + ', :) ' + ''.join(['{:16.7f}'.format(x) for x in jacobian[6*idx:6*(idx+1)]]))

print('')
