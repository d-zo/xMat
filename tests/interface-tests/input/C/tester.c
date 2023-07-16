#include <stdio.h>
#include <stdlib.h>


extern void Dot_Values_CInter(const char *materialname, const int *nparams,
   double materialparameters[], const int *nstatevar, double statevariables[],
   const int *ncomponents, double oldstress[], double oldstrain[], double timeincrement[],
   double totaltime[], double dotstress[], double dotstate[], double jacobian[]);

extern void xmat_console_CInter(const char *materialname, const int *nparams,
   double materialparameters[], const int *nstatevar, double statevariables[],
   const int *ncomponents, double oldstress[], double oldstrain[], double timeincrement[],
   double totaltime[], double newstress[], double newstate[], double jacobian[]);


typedef enum {
   num_components = 6,
   num_materialparameters = 16,
   num_statevariables = 20
} fixedLengths;


double *MemoryAllocationAndZeroArray(const int num_elements) {
   double *mem_newarray = calloc(num_elements, sizeof(double));
   return mem_newarray;
}


double *MemoryAllocationAndArrayAssignment(const int num_elements, double *originaldata) {
   double *mem_newarray = malloc(num_elements*sizeof(double));
   for (int idx = 0; idx < num_elements; idx++) {
      mem_newarray[idx] = originaldata[idx];
   }
   return mem_newarray;
}


int main() {
   char materialname[80] = "HYPO-VW96_Test      ";

   double p_val_materialparameters[num_materialparameters] = {0.5777, 0.0, 4000000.0, 0.27, 0.677,
      1.054, 1.212, 0.14, 2.5, 1.1, 2.2, 0.0001, 0.1, 5.5, 0.0, 0.8};
   double p_val_statevariables[num_statevariables] = {0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   double p_val_oldstress[num_components] = {-100.0, -250.0, -250.0, 50.0, 10.0, 0.0};
   double p_val_oldstrain[num_components] = {-0.00005, 0.000015, 0.000015, 0.0, 0.0, 0.0};

   double p_timeincrement = 0.001;
   double p_totaltime = 0.0;

   // Save array length as variables (so they can be referenced by an address)
   int r_num_materialparameters = num_materialparameters;
   int r_num_statevariables = num_statevariables;
   int r_num_components = num_components;

   // Allocate contiguous memory for input and output arrays
   double *p_materialparameters = MemoryAllocationAndArrayAssignment(num_materialparameters,
      p_val_materialparameters);
   double *p_statevariables = MemoryAllocationAndArrayAssignment(num_statevariables,
      p_val_statevariables);
   double *p_oldstress = MemoryAllocationAndArrayAssignment(num_components, p_val_oldstress);
   double *p_oldstrain = MemoryAllocationAndArrayAssignment(num_components, p_val_oldstrain);

   double *p_dotstress = MemoryAllocationAndZeroArray(num_components);
   double *p_dotstate = MemoryAllocationAndZeroArray(num_statevariables);
   double *p_jacobian = MemoryAllocationAndZeroArray(num_components*num_components);

   // Call the interface function and do the calculation
   Dot_Values_CInter(materialname, &r_num_materialparameters, p_materialparameters,
      &r_num_statevariables, p_statevariables, &r_num_components, p_oldstress, p_oldstrain,
      &p_timeincrement, &p_totaltime, p_dotstress, p_dotstate, p_jacobian);

   printf("dotstress = ");
   for (int idx = 0; idx < num_components; idx++) {
      printf("%16.7f", p_dotstress[idx]);
   }
   printf("\n");

   printf("dotstate = ");
   for (int idx = 0; idx < num_statevariables; idx++) {
      printf("%16.7f", p_dotstate[idx]);
   }
   printf("\n");

   printf("jacobian = ");
   for (int idx = 0; idx < num_components*num_components; idx++) {
      if ((idx % num_components) == 0) {
         printf("\n( %i, :) ", (int)(1+idx/num_components));
      }
      printf("%16.7f", p_jacobian[idx]);
   }
   printf("\n\n");

   double *p_newstress = p_dotstress;
   double *p_newstate = p_dotstate;

   xmat_console_CInter(materialname, &r_num_materialparameters, p_materialparameters,
      &r_num_statevariables, p_statevariables, &r_num_components, p_oldstress, p_oldstrain,
      &p_timeincrement, &p_totaltime, p_newstress, p_newstate, p_jacobian);

   printf("newstress = ");
   for (int idx = 0; idx < num_components; idx++) {
      printf("%16.7f", p_newstress[idx]);
   }
   printf("\n");

   printf("newstate = ");
   for (int idx = 0; idx < num_statevariables; idx++) {
      printf("%16.7f", p_newstate[idx]);
   }
   printf("\n");

   printf("jacobian = ");
   for (int idx = 0; idx < num_components*num_components; idx++) {
      if ((idx % num_components) == 0) {
         printf("\n( %i, :) ", (int)(1+idx/num_components));
      }
      printf("%16.7f", p_jacobian[idx]);
   }
   printf("\n\n");

   // Free allocated memory
   free(p_materialparameters);
   free(p_statevariables);
   free(p_oldstress);
   free(p_oldstrain);
   free(p_dotstress);
   free(p_dotstate);
   free(p_jacobian);

   return 0;
}
