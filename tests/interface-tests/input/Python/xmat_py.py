#!/usr/bin/python3
# -*- coding: utf-8 -*-


class xMat(object):
   """Wrapper class for Python to access the C interface of the xMat library.
   """
   # ----------------------------------------------------------------------------------------------
   def __init__(self, lib_path):
      import os
      #
      self.lib = None
      #
      loaded = False
      if (os.name == 'posix'):
         from ctypes import cdll
         #
         try:
            self.lib = cdll.LoadLibrary(lib_path)
            loaded = True
         except:
            print('Could not load <xmat>.so library (check path/name)')
      elif (os.name == 'nt'):
         from ctypes import WinDLL
         #
         try:
            self.lib = WinDLL(lib_path)
            loaded = True
         except:
            print('Could not load <xmat>.dll library (check path/name)')
      else:
         print('OS not supported')
      #
      if (loaded):
         from ctypes import c_char_p, c_int, c_double, POINTER
         self.lib.Dot_Values_CInter.restype = None
         self.lib.Dot_Values_CInter.argtypes = [  # --- Expected C Interface
            c_char_p,                             # const char *materialname,
            POINTER(c_int),                       # const int *nparams,
            POINTER(c_double),                    # double materialparameters[],
            POINTER(c_int),                       # const int *nstatevar,
            POINTER(c_double),                    # double statevariables[],
            POINTER(c_int),                       # const int *ncomponents,
            POINTER(c_double),                    # double oldstress[],
            POINTER(c_double),                    # double oldstrain[],
            POINTER(c_double),                    # double timeincrement[],
            POINTER(c_double),                    # double totaltime[],
            POINTER(c_double),                    # double dotstress[],
            POINTER(c_double),                    # double dotstate[],
            POINTER(c_double)                     # double jacobian[]
         ]
         self.lib.xmat_console_CInter.restype = None
         self.lib.xmat_console_CInter.argtypes = [# --- Expected C Interface
            c_char_p,                             # const char *materialname,
            POINTER(c_int),                       # const int *nparams,
            POINTER(c_double),                    # double materialparameters[],
            POINTER(c_int),                       # const int *nstatevar,
            POINTER(c_double),                    # double statevariables[],
            POINTER(c_int),                       # const int *ncomponents,
            POINTER(c_double),                    # double oldstress[],
            POINTER(c_double),                    # double oldstrain[],
            POINTER(c_double),                    # double timeincrement[],
            POINTER(c_double),                    # double totaltime[],
            POINTER(c_double),                    # double newstress[],
            POINTER(c_double),                    # double newstate[],
            POINTER(c_double)                     # double jacobian[]
         ]
   #
   # ----------------------------------------------------------------------------------------------
   def _PrepareInput(self, materialname, materialparameters, statevariables, oldstress,
      oldstrain, timeincrement, totaltime):
      import numpy
      from ctypes import c_double, c_char_p, c_int
      #
      c_materialname = c_char_p(materialname.encode())
      #
      nparams = len(materialparameters)
      c_nparams = c_int(nparams)
      cont_materialparameters = numpy.ascontiguousarray(materialparameters, dtype=numpy.double)
      c_materialparameters = (c_double * nparams)(*cont_materialparameters)
      #
      nstatevar = len(statevariables)
      c_nstatevar = c_int(nstatevar)
      cont_statevariables = numpy.ascontiguousarray(statevariables, dtype=numpy.double)
      c_statevariables = (c_double * nstatevar)(*cont_statevariables)
      #
      ncomponents = len(oldstress)
      c_ncomponents = c_int(ncomponents)
      cont_oldstress = numpy.ascontiguousarray(oldstress, dtype=numpy.double)
      c_oldstress = (c_double * ncomponents)(*cont_oldstress)
      cont_oldstrain = numpy.ascontiguousarray(oldstrain, dtype=numpy.double)
      c_oldstrain = (c_double * ncomponents)(*cont_oldstrain)
      #
      c_timeincrement = c_double(timeincrement)
      c_totaltime = c_double(totaltime)
      #
      # Preparation of return values
      outstress = [0.0 for x in range(ncomponents)]
      outstress = numpy.ascontiguousarray(outstress, dtype=numpy.double)
      double_ncomponents = c_double * ncomponents
      c_outstress = double_ncomponents(*outstress)
      #
      outstate = [0.0 for x in range(nstatevar)]
      outstate = numpy.ascontiguousarray(outstate, dtype=numpy.double)
      double_nstatevar = c_double * nstatevar
      c_outstate = double_nstatevar(*outstate)

      jacobian = [0.0 for x in range(ncomponents*ncomponents)]
      jacobian = numpy.ascontiguousarray(jacobian, dtype=numpy.double)
      double_n2components = c_double * (ncomponents*ncomponents)
      c_jacobian = double_n2components(*jacobian)
      #
      return [c_materialname, c_nparams, c_materialparameters,
         c_nstatevar, c_statevariables, c_ncomponents, c_oldstress, c_oldstrain,
         c_timeincrement, c_totaltime, c_outstress, c_outstate, c_jacobian]
   #
   # ----------------------------------------------------------------------------------------------
   def Dot_Values_CInter(self, materialname, materialparameters, statevariables, oldstress,
      oldstrain, timeincrement, totaltime):
      from ctypes import byref, c_double, POINTER, cast
      #
      if (self.lib is None):
         print('No library loaded')
         return [None, None, None]
      #
      c_materialname, c_nparams, c_materialparameters, c_nstatevar, c_statevariables, \
      c_ncomponents, c_oldstress, c_oldstrain, c_timeincrement, c_totaltime, c_dotstress, \
      c_dotstate, c_jacobian = self._PrepareInput(materialname=materialname,
         materialparameters=materialparameters, statevariables=statevariables, oldstress=oldstress,
         oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime)
      #
      # Call the routine
      self.lib.Dot_Values_CInter(c_materialname, byref(c_nparams), c_materialparameters,
         byref(c_nstatevar), c_statevariables, byref(c_ncomponents), c_oldstress, c_oldstrain,
         byref(c_timeincrement), byref(c_totaltime), c_dotstress, c_dotstate, c_jacobian)
      #
      nparams = len(materialparameters)
      nstatevar = len(statevariables)
      ncomponents = len(oldstress)
      #
      # Assign return values to python lists
      temp_dotstress = cast(c_dotstress, POINTER(c_double*ncomponents))
      dotstress = [temp_dotstress.contents[idx] for idx in range(ncomponents)]
      #
      temp_dotstate = cast(c_dotstate, POINTER(c_double*nstatevar))
      dotstate = [temp_dotstate.contents[idx] for idx in range(nstatevar)]
      #
      temp_jacobian = cast(c_jacobian, POINTER(c_double*(ncomponents*ncomponents)))
      jacobian = [temp_jacobian.contents[idx] for idx in range(ncomponents*ncomponents)]
      #
      return [dotstress, dotstate, jacobian]
   #
   # ----------------------------------------------------------------------------------------------
   def xmat_console_CInter(self, materialname, materialparameters, statevariables, oldstress,
      oldstrain, timeincrement, totaltime):
      from ctypes import byref, c_double, POINTER, cast
      #
      if (self.lib is None):
         print('No library loaded')
         return [None, None, None]
      #
      c_materialname, c_nparams, c_materialparameters, c_nstatevar, c_statevariables, \
      c_ncomponents, c_oldstress, c_oldstrain, c_timeincrement, c_totaltime, c_newstress, \
      c_newstate, c_jacobian = self._PrepareInput(materialname=materialname,
         materialparameters=materialparameters, statevariables=statevariables, oldstress=oldstress,
         oldstrain=oldstrain, timeincrement=timeincrement, totaltime=totaltime)
      #
      # Call the routine
      self.lib.xmat_console_CInter(c_materialname, byref(c_nparams), c_materialparameters,
         byref(c_nstatevar), c_statevariables, byref(c_ncomponents), c_oldstress, c_oldstrain,
         byref(c_timeincrement), byref(c_totaltime), c_newstress, c_newstate, c_jacobian)
      #
      nparams = len(materialparameters)
      nstatevar = len(statevariables)
      ncomponents = len(oldstress)
      #
      # Assign return values to python lists
      temp_newstress = cast(c_newstress, POINTER(c_double*ncomponents))
      newstress = [temp_newstress.contents[idx] for idx in range(ncomponents)]
      #
      temp_newstate = cast(c_newstate, POINTER(c_double*nstatevar))
      newstate = [temp_newstate.contents[idx] for idx in range(nstatevar)]
      #
      temp_jacobian = cast(c_jacobian, POINTER(c_double*(ncomponents*ncomponents)))
      jacobian = [temp_jacobian.contents[idx] for idx in range(ncomponents*ncomponents)]
      #
      return [newstress, newstate, jacobian]
#
