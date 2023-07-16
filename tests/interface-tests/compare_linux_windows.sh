#!/bin/bash

echo "Comparison remark: No diff output means identical except trailing carriage return"

echo -e "\nComparing Windows/Linux"
for testrun in "C" "Fortran" "Fortran-UMAT" "Fortran-USER_MOD" "Fortran-VUMAT" "Python"; do
   echo " - Comparing output/${testrun}/output.lnx output/${testrun}/output.win"
   diff --strip-trailing-cr output/${testrun}/output.lnx output/${testrun}/output.win
done

echo -e "\nComparing to C reference"
for testrun in "Fortran" "Python"; do
   echo " - Comparing output/C/output.lnx output/${testrun}/output.lnx"
   diff --strip-trailing-cr output/C/output.lnx output/${testrun}/output.lnx
done

echo -e "\nComparing to last part of C reference"
for testrun in "Fortran-UMAT" "Fortran-USER_MOD"; do
   echo " - Comparing (part of) output/C/output.lnx output/${testrun}/output.lnx"
   diff --strip-trailing-cr <(tail -n +11 output/C/output.lnx) output/${testrun}/output.lnx
done
for testrun in "Matlab" "Octave"; do
   echo " - Comparing (part of) output/C/output.lnx output/${testrun}/output"
   diff --strip-trailing-cr <(tail -n +11 output/C/output.lnx) output/${testrun}/output
done

echo -e "\nComparing to part of C reference"
echo " - Comparing (part of) output/C/output.lnx output/Fortran-VUMAT/output.lnx"
diff --strip-trailing-cr <(tail -n +11 output/C/output.lnx | head -n 2) output/Fortran-VUMAT/output.lnx

echo -e "\nComparing to first part of C reference"
echo " - Comparing (part of) output/C/output.lnx output/Simulink/output"
diff --strip-trailing-cr <(head -n 12 output/C/output.lnx) output/Simulink/output

exit 0
