#!/bin/bash

cd src/Math_Operations/;
for srcfile in common/*; do
   rmfile="$(basename ${srcfile})";
   for dir in "mat66" "mat99" "tens3333"; do
      if [ -e ${dir}/${rmfile} ]; then
         rm ${dir}/${rmfile};
      fi
      ln -s "../${srcfile}" ${dir}/${rmfile};
   done
done

exit 0;
