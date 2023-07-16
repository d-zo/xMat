% Run this script in MATLAB after compiling xmat.f and xmat_sl.c as a C MEX S-function
% and creating the modell xmat_test.slx

% Start test system in Simulink and run it
open_system('xmat_test.slx');
sim_res = sim('xmat_test', 'SimulationMode', 'normal', 'CaptureErrors', 'on')
disp(sim_res.ErrorMessage);

disp('-- Done --');

dotstate = sim_res.dotstate(end, :);
dotstress = sim_res.dotstress(end, :);
jacobian = sim_res.jacobian(end, :);
newstate = sim_res.newstate(end, :);
newstress = sim_res.newstress(end, :);

fid = fopen('output', 'w');
fprintf(fid, 'dotstress = ');
fprintf(fid, '%16.7f', dotstress);
fprintf(fid, '\n');

fprintf(fid, 'dotstate = ');
fprintf(fid, '%16.7f', dotstate);
fprintf(fid, '\n');

fprintf(fid, 'jacobian = \n');
for idx = 1:6
   fprintf(fid, '( %i, :) ', idx);
   fprintf(fid, '%16.7f', jacobian(((idx-1)*6+1):6*idx));
   fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, 'newstress = ');
fprintf(fid, '%16.7f', newstress);
fprintf(fid, '\n');

fprintf(fid, 'newstate = ');
fprintf(fid, '%16.7f', newstate);
fprintf(fid, '\n');
fclose(fid);
