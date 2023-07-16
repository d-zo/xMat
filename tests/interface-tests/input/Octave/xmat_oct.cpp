#include <octave/oct.h>


extern "C" {
   void xmat_console_CInter(const char *materialname, const int *nparams,
      double materialparameters[], const int *nstatevar, double statevariables[],
      const int *ncomponents, double oldstress[], double oldstrain[], double timeincrement[],
      double totaltime[], double newstress[], double newstate[], double jacobian[]);
}


DEFUN_DLD(xmat_oct, args, , "Calling xMat from Octave") {
   octave_value_list retval;
   if (args.length () != 7) {
      print_usage();
      return retval;
   }

   std::string s_materialname = args(0).string_value();
   const char* materialname = s_materialname.c_str();

   NDArray array_materialparameters = args(1).array_value();
   octave_idx_type o_num_materialparameters = array_materialparameters.numel();
   int r_num_materialparameters = o_num_materialparameters;
   double *p_materialparameters = array_materialparameters.fortran_vec();

   NDArray array_statevariables = args(2).array_value();
   octave_idx_type o_num_statevariables = array_statevariables.numel();
   int num_statevariables = o_num_statevariables;
   double *p_statevariables = array_statevariables.fortran_vec();

   NDArray array_oldstress = args(3).array_value();
   octave_idx_type o_num_components = array_oldstress.numel();
   int num_components = o_num_components;
   double *p_oldstress = array_oldstress.fortran_vec();

   NDArray array_oldstrain = args(4).array_value();
   double *p_oldstrain = array_oldstrain.fortran_vec();

   octave_value val_timeincrement = args(5);
   double p_timeincrement = val_timeincrement.double_value();

   octave_value val_totaltime = args(6);
   double p_totaltime = val_totaltime.double_value();

   OCTAVE_LOCAL_BUFFER(double, p_newstress, num_components);
   OCTAVE_LOCAL_BUFFER(double, p_newstate, num_statevariables);
   OCTAVE_LOCAL_BUFFER(double, p_jacobian, num_components*num_components);

   // Call the interface function and do the calculation
   xmat_console_CInter(materialname, &r_num_materialparameters, p_materialparameters,
      &num_statevariables, p_statevariables, &num_components, p_oldstress, p_oldstrain,
      &p_timeincrement, &p_totaltime, p_newstress, p_newstate, p_jacobian);

   dim_vector dim_stress(num_components, 1);
   NDArray newstress(dim_stress);
   for (int idx = 0; idx < num_components; idx++) {
      newstress(idx) = p_newstress[idx];
   }
   retval.append(newstress);

   dim_vector dim_states(num_statevariables, 1);
   NDArray newstate(dim_states);
   for (int idx = 0; idx < num_statevariables; idx++) {
      newstate(idx) = p_newstate[idx];
   }
   retval.append(newstate);

   dim_vector dim_jacobi(num_components, num_components);
   NDArray jacobian(dim_jacobi);
   for (int idx = 0; idx < num_components; idx++) {
      for (int jdx = 0; jdx < num_components; jdx++) {
         // Transpose matrix while assigning components back
         jacobian(jdx, idx) = p_jacobian[idx*num_components+jdx];
      }
   }
   retval.append(jacobian);

   return retval;
}
