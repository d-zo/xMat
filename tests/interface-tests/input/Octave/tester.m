% Run this script in GNU Octave after compiling xmat.f as an OCT file

materialname = 'HYPO-VW96_Test';
materialparameters = [0.5777, 0.0, 4000000.0, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5, ...
                      1.1, 2.2, 0.0001, 0.1, 5.5, 0.0, 0.8];
statevariables = [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
oldstress = [-100.0, -250.0, -250.0, 50.0, 10.0, 0.0];
oldstrain = [-0.00005, 0.000015, 0.000015, 0.0, 0.0, 0.0];
timeincrement = 0.001;
totaltime = 0.0;

start=tic();
[newstress, newstate, jacobian] = xmat_oct(materialname, materialparameters, statevariables, ...
                                           oldstress, oldstrain, timeincrement, totaltime);
toc(start)

fid = fopen('output', 'w');
fprintf(fid, 'newstress = ');
fprintf(fid, '%16.7f', newstress);
fprintf(fid, '\n');

fprintf(fid, 'newstate = ');
fprintf(fid, '%16.7f', newstate);
fprintf(fid, '\n');

fprintf(fid, 'jacobian = \n');
for idx = 1:6
   fprintf(fid, '( %i, :) ', idx);
   fprintf(fid, '%16.7f', jacobian(((idx-1)*6+1):6*idx));
   fprintf(fid, '\n');
end
fprintf(fid, '\n');
fclose(fid);
