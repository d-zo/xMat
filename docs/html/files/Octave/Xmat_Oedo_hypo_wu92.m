% Run this script in GNU Octave after compiling xmat.f as an OCT file

maxiter = 50000;

materialname = 'Hypo-Wu92_Karlsruher Sand dicht 01';
voidratio = 0.65;
materialparameters = [-106.5, -801.5, -797.1, 1077.7];
statevariables = [voidratio, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

oedo_pressure = [-10.0, -100.0, -50.0, -200.0];
timeincrement = 0.0001;
totaltime = 0.0001;

stress = [oedo_pressure(1)/2, oedo_pressure(1), oedo_pressure(1)/2, 0.0, 0.0, 0.0];
strain = [0.0, -0.05, 0.0, 0.0, 0.0, 0.0]*timeincrement;

oldstress = stress;
oldstate = statevariables;

idx = 1;
breakall = false;
numbersign = 1.0;

fid = fopen('Xmat_Oedo_Octave.csv', 'w');

fprintf(fid, '%13.6f   %13.6f\n', -oldstress(2), voidratio);

start=tic();
for istep = 2:numel(oedo_pressure)
   numbersign = (-1)^istep;
   oldstrain = numbersign*strain;
   while (true)
      if (numbersign*oldstress(2) < numbersign*oedo_pressure(istep))
         break;
      end

      [newstress, newstate, jacobian] = xmat_oct(materialname, materialparameters, oldstate, ...
         oldstress, oldstrain, timeincrement, totaltime);
   
      voidratio = newstate(1);
      oldstress = newstress;
      oldstate = newstate;
      
      fprintf(fid, '%13.6f   %13.6f\n', -oldstress(2), voidratio);

      idx += 1;
      if (idx > maxiter)
         breakall = true;
         break;
      end
   end
   if (breakall)
      break;
   end
end
printf("Xmat_Oedo_Octave: %6.3fs\n", toc(start));

fclose(fid);
