#define S_FUNCTION_NAME  xmat_sl
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO


extern void Dot_Values_CInter(const char *materialname, const int *nparams,
   double materialparameters[], const int *nstatevar, double statevariables[],
   const int *ncomponents, double oldstress[], double oldstrain[], double timeincrement[],
   double totaltime[], double dotstress[], double dotstate[], double jacobian[]);


typedef enum {
   idx_dotstress = 0,
   idx_dotstate,
   idx_jacobian,
   num_outputparams
} outputParam;


typedef enum {
   idx_materialname = 0,
   idx_materialparameters,
   idx_statevariables,
   idx_oldstress,
   idx_oldstrain,
   idx_timeincrement,
   idx_totaltime,
   num_inputparams
} inputParam;


// Currently fixed port widths are used
typedef enum {
   num_components = 6,
   num_materialparameters = 16,
   num_statevariables = 20
} fixedPortWidths;


static void mdlInitializeSizes(SimStruct *S) {
   // Set the number of input ports
   if (!ssSetNumInputPorts(S, num_inputparams)) {
      return;
   }
   // Set the number of output ports
   if (!ssSetNumOutputPorts(S, num_outputparams)) {
      return;
   }
   // Number of expected parameters
   ssSetNumSFcnParams(S, 0);

   int in_dataportwidth[num_inputparams] = {DYNAMICALLY_SIZED, num_materialparameters,
      num_statevariables, num_components, num_components, 1, 1};
   int out_dataportwidth[num_inputparams] = {num_components, num_statevariables,
      num_components*num_components};

   // Set the number of work vectors (contiguous allocated memory provided for Fortran output values)
   ssSetNumDWork(S, num_outputparams);
   for (int idx = 0; idx < num_outputparams; idx++) {
      ssSetDWorkWidth(S, idx, out_dataportwidth[idx]);
      ssSetDWorkDataType(S, idx, SS_DOUBLE);
   }

   // Prepare all input ports
   for (int idx = 0; idx < num_inputparams; idx++) {
      if (idx == idx_materialname) {
         ssSetInputPortDataType(S, idx, SS_UINT8);
      }
      else {
         ssSetInputPortDataType(S, idx, SS_DOUBLE);
      }
      ssSetInputPortWidth(S, idx, in_dataportwidth[idx]);
      ssSetInputPortComplexSignal(S, idx, COMPLEX_NO);
      ssSetInputPortDirectFeedThrough(S, idx, 1);
      ssSetInputPortAcceptExprInRTW(S, idx, 0);
      ssSetInputPortOverWritable(S, idx, 0);
      ssSetInputPortOptimOpts(S, idx, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortRequiredContiguous(S, idx, 1);
   }

   // Prepare all output ports
   for (int idx = 0; idx < num_outputparams; idx++) {
      ssSetOutputPortDataType(S, idx, SS_DOUBLE);
      ssSetOutputPortWidth(S, idx, out_dataportwidth[idx]);
      ssSetOutputPortComplexSignal(S, idx, COMPLEX_NO);
      ssSetOutputPortOptimOpts(S, idx, SS_REUSABLE_AND_LOCAL);
      ssSetOutputPortOutputExprInRTW(S, idx, 0);
   }

   // Register reserved identifiers to avoid name conflict
   if (ssRTWGenIsCodeGen(S) || ssGetSimMode(S)==SS_SIMMODE_EXTERNAL) {
      ssRegMdlInfo(S, "dot_values", MDL_INFO_ID_RESERVED, 0, 0, ssGetPath(S));
   }  

   // This S-function can be used in referenced model simulating in normal mode
   ssSetModelReferenceNormalModeSupport(S, MDL_START_AND_MDL_PROCESS_PARAMS_OK);

   // Set the number of sample time
   ssSetNumSampleTimes(S, 1);

   // Set the compliance for the operating point save/restore.
   ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

   ssSetArrayLayoutForCodeGen(S, SS_ALL);

   // Set the Simulink version this S-Function has been generated in
   ssSetSimulinkVersionGeneratedIn(S, "10.0");

   // All options are documented in matlabroot/simulink/include/simstruc.h
   ssSetOptions(S,
      SS_OPTION_USE_TLC_WITH_ACCELERATOR |
      SS_OPTION_CAN_BE_CALLED_CONDITIONALLY |
      SS_OPTION_EXCEPTION_FREE_CODE |
      SS_OPTION_WORKS_WITH_CODE_REUSE |
      SS_OPTION_SFUNCTION_INLINED_FOR_RTW |
      SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME
   );
}


static void mdlInitializeSampleTimes(SimStruct *S) {
   ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
   ssSetOffsetTime(S, 0, FIXED_IN_MINOR_STEP_OFFSET);

   #if defined(ssSetModelReferenceSampleTimeDefaultInheritance)
   ssSetModelReferenceSampleTimeDefaultInheritance(S);
   #endif
}


static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T portIndex, const DimsInfo_T *dimsInfo) {
   // Set input port dimension
   if (!ssSetInputPortDimensionInfo(S, portIndex, dimsInfo)) { 
      return;
   }
}


static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T portIndex, const DimsInfo_T *dimsInfo) {
   // Set output port dimension
   if(!ssSetOutputPortDimensionInfo(S, portIndex, dimsInfo)) {
      return;
   }
}


static void mdlSetDefaultPortDimensionInfo(SimStruct *S) {
   // Default dimension for all input ports
   for (int idx = 1; idx <= num_inputparams; idx++) {
      if (ssGetInputPortWidth(S, idx) == DYNAMICALLY_SIZED) {
         ssSetInputPortWidth(S, idx, 1);
      }
   }
}


static void mdlStart(SimStruct *S) {
   // Allocate the memory availabe for DWork data and initialize it with zeros
   double* mem_dotstress = (double*) ssGetDWork(S, 0);
   for (int idx = 0; idx < num_components; idx++) {
      mem_dotstress[idx] = 0.0;
   }
   double* mem_dotstate = (double*) ssGetDWork(S, 1);
   for (int idx = 0; idx < num_statevariables; idx++) {
      mem_dotstate[idx] = 0.0;
   }
   double* mem_jacobian = (double*) ssGetDWork(S, 2);
   for (int idx = 0; idx < num_components*num_components; idx++) {
      mem_jacobian[idx] = 0.0;
   }
}


static void mdlOutputs(SimStruct *S, int_T tid) {
   // Get access to input, output and DWork data
   char* f_materialname = (char*) ssGetInputPortSignal(S, idx_materialname);
   double* p_materialparameters = (double*) ssGetInputPortSignal(S, idx_materialparameters);
   double* p_statevariables = (double*) ssGetInputPortSignal(S, idx_statevariables);
   double* p_oldstress = (double*) ssGetInputPortSignal(S, idx_oldstress);
   double* p_oldstrain = (double*) ssGetInputPortSignal(S, idx_oldstrain);
   double* p_timeincrement = (double*) ssGetInputPortSignal(S, idx_timeincrement);
   double* p_totaltime = (double*) ssGetInputPortSignal(S, idx_totaltime);

   double* p_dotstress = (double*) ssGetOutputPortSignal(S, idx_dotstress);
   double* p_dotstate = (double*) ssGetOutputPortSignal(S, idx_dotstate);
   double* p_jacobian = (double*) ssGetOutputPortSignal(S, idx_jacobian);

   double* mem_dotstress = (double*) ssGetDWork(S, 0);
   double* mem_dotstate = (double*) ssGetDWork(S, 1);
   double* mem_jacobian = (double*) ssGetDWork(S, 2);

   // Save array length as variables (so they have an address as a reference)
   int r_num_materialparameters = num_materialparameters;
   int r_num_statevariables = num_statevariables;
   int r_num_components = num_components;

   // Call the xMat interface function
   Dot_Values_CInter(f_materialname, &r_num_materialparameters, p_materialparameters,
      &r_num_statevariables, p_statevariables, &r_num_components, p_oldstress, p_oldstrain,
      p_timeincrement, p_totaltime, mem_dotstress, mem_dotstate, mem_jacobian);

   // Reassign output values
   for (int idx = 0; idx < num_components; idx++) {
      p_dotstress[idx] = mem_dotstress[idx];
   }
   for (int idx = 0; idx < num_statevariables; idx++) {
      p_dotstate[idx] = mem_dotstate[idx];
   }
   for (int idx = 0; idx < num_components*num_components; idx++) {
      p_jacobian[idx] = mem_jacobian[idx];
   }
}


static void mdlTerminate(SimStruct *S) {
   // This method must exist, but nothing has to be done for this block
}


#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
