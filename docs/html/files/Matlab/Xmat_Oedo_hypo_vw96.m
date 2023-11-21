% Run this script in MATLAB after compiling xmat.f as an MEX file

maxiter = 50000;

materialname = 'Hypo-VW96_Hostun Sand 01';
voidratio = 0.7;
nu = 0.25;
materialparameters = [0.576, nu, 1.0e6, 0.29, 0.63, 1.0, 1.15, 0.13, 2.0, ...
                      2.0, 5.0, 1.0e-4, 0.5, 6.0, 0.0, 0.0];
intergranular_strain = [0.0, -0.0001, 0.0, 0.0, 0.0, 0.0];
statevariables = [voidratio, intergranular_strain(1), intergranular_strain(2), intergranular_strain(3), ...
                  intergranular_strain(4), intergranular_strain(4), intergranular_strain(5), ...
                  intergranular_strain(5), intergranular_strain(6), intergranular_strain(6), ...
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

fid = fopen('Xmat_Oedo_Matlab.csv', 'w');

fprintf(fid, '%13.6f   %13.6f\n', -oldstress(2), voidratio);

start=tic();
for istep = 2:numel(oedo_pressure)
   numbersign = (-1)^istep;
   oldstrain = numbersign*strain;
   while (true)
      if (numbersign*oldstress(2) < numbersign*oedo_pressure(istep))
         break;
      end

      [newstress, newstate, jacobian] = xmat(materialname, materialparameters, oldstate, ...
         oldstress, oldstrain, timeincrement, totaltime);
   
      voidratio = newstate(1);
      oldstress = newstress;
      oldstate = newstate;
      
      fprintf(fid, '%13.6f   %13.6f\n', -oldstress(2), voidratio);

      idx = idx+1;
      if (idx > maxiter)
         breakall = true;
         break;
      end
   end
   if (breakall)
      break;
   end
end
sprintf("Xmat_Oedo_Matlab: %6.3fs\n", toc(start));

fclose(fid);
