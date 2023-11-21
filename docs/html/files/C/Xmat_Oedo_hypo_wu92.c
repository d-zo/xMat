#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>


extern void xmat_console_CInter(const char *materialname, const int *nparams,
   double materialparameters[], const int *nstatevar, double statevariables[],
   const int *ncomponents, double oldstress[], double oldstrain[], double timeincrement[],
   double totaltime[], double newstress[], double newstate[], double jacobian[]);


typedef enum {
   num_components = 6,
   num_materialparameters = 4,
   num_statevariables = 20,
   num_oedo_pressure = 4
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
   char materialname[80] = "Hypo-Wu92_Karlsruher Sand dicht 01";

   double voidratio = 0.65;
   double p_val_materialparameters[num_materialparameters] = {-106.5, -801.5, -797.1, 1077.7};
   double p_val_statevariables[num_statevariables] = {voidratio, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   double oedo_pressure[num_oedo_pressure] = {-10.0, -100.0, -50.0, -200.0};

   double p_timeincrement = 0.0001;
   double p_totaltime = 0.0001;

   double p_val_stress[num_components] = {oedo_pressure[0]/2, oedo_pressure[0], oedo_pressure[0]/2, 0.0, 0.0, 0.0};
   double p_val_strain[num_components] = {0.0, -0.05*p_timeincrement, 0.0, 0.0, 0.0, 0.0};

   const int maxiter = 50000;
   int idx = 1;
   bool breakall = false;
   double numbersign = 1.0;

   // Save array length as variables (so they can be referenced by an address)
   int r_num_materialparameters = num_materialparameters;
   int r_num_statevariables = num_statevariables;
   int r_num_components = num_components;


   // Allocate contiguous memory for input and output arrays
   double *p_materialparameters = MemoryAllocationAndArrayAssignment(num_materialparameters,
      p_val_materialparameters);
   double *p_statevariables = MemoryAllocationAndArrayAssignment(num_statevariables,
      p_val_statevariables);
   double *p_oldstress = MemoryAllocationAndArrayAssignment(num_components, p_val_stress);
   double *p_oldstrain = MemoryAllocationAndArrayAssignment(num_components, p_val_strain);

   double *p_newstress = MemoryAllocationAndZeroArray(num_components);
   double *p_newstate = MemoryAllocationAndZeroArray(num_statevariables);
   double *p_jacobian = MemoryAllocationAndZeroArray(num_components*num_components);

   clock_t start_time, end_time; 
   start_time = clock(); 

   FILE *file_handle = fopen("Xmat_Oedo_C.csv", "w");
   fprintf(file_handle, "%13.6f   %13.6f\n", -p_oldstress[1], voidratio);

   for (int istep = 1; istep < num_oedo_pressure; istep++) {
      numbersign = pow(-1, istep+1);
      for (int icomp = 0; icomp < num_components; icomp++) {
         p_oldstrain[icomp] = numbersign*p_val_strain[icomp];
      }
      while (true) {
         if (numbersign*p_oldstress[1] < numbersign*oedo_pressure[istep]) {
            break;
         }

         xmat_console_CInter(materialname, &r_num_materialparameters, p_materialparameters,
            &r_num_statevariables, p_statevariables, &r_num_components, p_oldstress, p_oldstrain,
            &p_timeincrement, &p_totaltime, p_newstress, p_newstate, p_jacobian);

         voidratio = p_newstate[0];
         for (int icomp = 0; icomp < num_components; icomp++) {
            p_oldstress[icomp] = p_newstress[icomp];
         }
         for (int istatev = 0; istatev < num_statevariables; istatev++) {
            p_statevariables[istatev] = p_newstate[istatev];
         }

         fprintf(file_handle, "%13.6f   %13.6f\n", -p_oldstress[1], voidratio);

         idx++;
         if (idx > maxiter) {
            breakall = true;
            break;
         }
      }
      if (breakall) {
         break;
      }
   }
   end_time = clock(); 
   double elapsed_time = ((double)(end_time - start_time))/CLOCKS_PER_SEC;
   printf("Xmat_Oedo_C: %6.3fs\n", elapsed_time);

   fclose(file_handle);

   // Free allocated memory
   free(p_materialparameters);
   free(p_statevariables);
   free(p_oldstress);
   free(p_oldstrain);
   free(p_newstress);
   free(p_newstate);
   free(p_jacobian);

   return 0;
}
