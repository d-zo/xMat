%% Create a new Simulink system xmat_test
% If any system named xmat_test is open, close it without saving any changes.
% Also remove xmat_test model file if it exists
if bdIsLoaded('xmat_test')
   bdclose('xmat_test')
end
if exist('xmat_test.slx', 'file')
   delete('xmat_test.slx')
end

materialname = '"HYPO-VW96_Test"';
materialparameters = [0.5777, 0.0, 4000000.0, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5, ...
                      1.1, 2.2, 0.0001, 0.1, 5.5, 0.0, 0.8];
statevariables = [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
stress = [-100.0, -250.0, -250.0, 50.0, 10.0, 0.0];
strain = [-0.00005, 0.000015, 0.000015, 0.0, 0.0, 0.0];
dt = 0.001;
totaltime =  0.0;

new_system('xmat_test');

% Input blocks
add_block('built-in/StringConstant', 'xmat_test/materialname', 'String', materialname, ...
   'Position', [205, 120, 320, 150]);
add_block('built-in/StringToASCII', 'xmat_test/String to ASCII', 'Position', [445, 120, 540, 150]);
add_block('built-in/Constant', 'xmat_test/materialparameters', 'Value', mat2str(materialparameters), ...
   'Position', [225, 180, 320, 210]);
add_block('built-in/Constant', 'xmat_test/statevariables', 'Value', mat2str(statevariables), ...
   'Position', [225, 240, 320, 270]);
add_block('built-in/Constant', 'xmat_test/stress', 'Value', mat2str(stress), ...
   'Position', [225, 300, 320, 330]);
add_block('built-in/Constant', 'xmat_test/strain', 'Value', mat2str(strain), ...
   'Position', [225, 360, 320, 390]);
add_block('built-in/Constant', 'xmat_test/dt', 'Value', num2str(dt), 'Position', [290, 420, 320, 450]);
add_block('built-in/Constant', 'xmat_test/totaltime', 'Value', num2str(totaltime), ...
   'Position', [290, 480, 320, 510]);

% xMat S-Function block
add_block('built-in/S-Function', 'xmat_test/S-Function of xMat', 'FunctionName', 'xmat_sl', ...
   'Position', [445, 220, 575, 410], 'BackgroundColor', '[0.569, 0.843, 0.263]');

% Simple one-step explicit Euler integration
add_block('hdlsllib/HDL Operations/Multiply-Add', 'xmat_test/Mult-Add-stress', ...
   'Position', [515, 0, 575, 50]);
add_block('hdlsllib/HDL Operations/Multiply-Add', 'xmat_test/Mult-Add-statev', ...
   'Position', [515, 580, 575, 630]);

% Output to workspace blocks
add_block('built-in/ToWorkspace', 'xmat_test/out_dotstress', 'VariableName', 'dotstress', ...
   'Position', [700, 140, 790, 170]);
add_block('built-in/ToWorkspace', 'xmat_test/out_newstress', 'VariableName', 'newstress', ...
   'Position', [700, 220, 790, 250]);
add_block('built-in/ToWorkspace', 'xmat_test/out_dotstate', 'VariableName', 'dotstate', ...
   'Position', [700, 300, 790, 330]);
add_block('built-in/ToWorkspace', 'xmat_test/out_newstate', 'VariableName', 'newstate', ...
   'Position', [700, 380, 790, 410]);
add_block('built-in/ToWorkspace', 'xmat_test/out_jacobian', 'VariableName', 'jacobian', ...
   'Position', [700, 460, 790, 490]);

% Display blocks
add_block('built-in/Display', 'xmat_test/new_stress', 'Position', [910, -100, 1000, 30]);
add_block('built-in/Display', 'xmat_test/dot_stress', 'Position', [910, 60, 1000, 190]);
add_block('built-in/Display', 'xmat_test/dot_state', 'Position', [910, 220, 1000, 350]);
add_block('built-in/Display', 'xmat_test/jacobian', 'Position', [910, 380, 1000, 510]);
add_block('built-in/Display', 'xmat_test/new_state', 'Position', [910, 540, 1000, 670]);

% Connections into xMat block
add_line('xmat_test', ['materialname', '/1'], ['String to ASCII', '/1']);
lh = add_line('xmat_test', ['String to ASCII', '/1'], ['S-Function of xMat', '/1']);
set(lh, 'Points', [540, 135; 560, 135; 560, 185; 430, 185; 430, 240; 445, 240]);
lh = add_line('xmat_test', ['materialparameters', '/1'], ['S-Function of xMat', '/2']);
set(lh, 'Points', [325, 195; 385, 195; 385, 265; 445, 265]);
lh = add_line('xmat_test', ['statevariables', '/1'], ['S-Function of xMat', '/3']);
set(lh, 'Points', [325, 255; 370, 255; 370, 290; 445, 290]);
add_line('xmat_test', ['stress', '/1'], ['S-Function of xMat', '/4']);
lh = add_line('xmat_test', ['strain', '/1'], ['S-Function of xMat', '/5']);
set(lh, 'Points', [325, 375; 385, 375; 385, 340; 445, 340]);
lh = add_line('xmat_test', ['dt', '/1'], ['S-Function of xMat', '/6']);
set(lh, 'Points', [325, 435; 400, 435; 400, 365; 445, 365]);
lh = add_line('xmat_test', ['totaltime', '/1'], ['S-Function of xMat', '/7']);
set(lh, 'Points', [325, 495; 430, 495; 430, 390; 445, 390]);
   
% Connections out of xMat block
lh = add_line('xmat_test', ['S-Function of xMat', '/1'], ['out_dotstress', '/1']);
set(lh, 'Points', [580, 250; 620, 250; 620, 155; 695, 155]);
lh = add_line('xmat_test', ['S-Function of xMat', '/1'], ['dot_stress', '/1']);
set(lh, 'Points', [580, 250; 620, 250; 620, 155; 680, 155; 680, 125; 905, 125]);
add_line('xmat_test', ['S-Function of xMat', '/2'], ['out_dotstate', '/1']);
lh = add_line('xmat_test', ['S-Function of xMat', '/2'], ['dot_state', '/1']);
set(lh, 'Points', [580, 315; 680, 315; 680, 285; 905, 285]);
lh = add_line('xmat_test', ['S-Function of xMat', '/3'], ['out_jacobian', '/1']);
set(lh, 'Points', [580, 380; 620, 380; 620, 475; 695, 475]);
lh = add_line('xmat_test', ['S-Function of xMat', '/3'], ['jacobian', '/1']);
set(lh, 'Points', [580, 380; 620, 380; 620, 475; 680, 475; 680, 445; 905, 445]);


% Connections to Multiply-Adds
lh = add_line('xmat_test', ['dt', '/1'], ['Mult-Add-stress', '/1']);
set(lh, 'Points', [325, 435; 400, 435; 400, 365; 400, 10; 510, 10]);
lh = add_line('xmat_test', ['S-Function of xMat', '/1'], ['Mult-Add-stress', '/2']);
set(lh, 'Points', [580, 250; 620, 250; 620, 155; 620, 75; 490, 75; 490, 25; 510, 25]);
lh = add_line('xmat_test', ['stress', '/1'], ['Mult-Add-stress', '/3']);
set(lh, 'Points', [325, 315; 415, 315; 415, 40; 510, 40]);
lh = add_line('xmat_test', ['Mult-Add-stress', '/1'], ['out_newstress', '/1']);
set(lh, 'Points', [575, 25; 640, 25; 640, 235; 695, 235]);
lh = add_line('xmat_test', ['Mult-Add-stress', '/1'], ['new_stress', '/1']);
set(lh, 'Points', [575, 25; 640, 25; 640, -35; 905, -35]);

lh = add_line('xmat_test', ['dt', '/1'], ['Mult-Add-statev', '/1']);
set(lh, 'Points', [325, 435; 400, 435; 400, 590; 510, 590]);
lh = add_line('xmat_test', ['S-Function of xMat', '/2'], ['Mult-Add-statev', '/2']);
set(lh, 'Points', [580, 315; 600, 315; 600, 555; 490, 555; 490, 605; 510, 605]);
lh = add_line('xmat_test', ['statevariables', '/1'], ['Mult-Add-statev', '/3']);
set(lh, 'Points', [325, 255; 370, 255; 370, 290; 370, 620; 510, 620]);
lh = add_line('xmat_test', ['Mult-Add-statev', '/1'], ['out_newstate', '/1']);
set(lh, 'Points', [575, 605; 640, 605; 640, 395; 695, 395]);
add_line('xmat_test', ['Mult-Add-statev', '/1'], ['new_state', '/1']);

% Annotations
note = Simulink.Annotation('xmat_test', 'Input_values');
set(note, 'FixedHeight', 'on', 'FixedWidth', 'on', 'Position', [175, 85, 355, 535], ...
   'ForegroundColor', '[0.45, 0.45, 0.45]', 'BackgroundColor', '[0.95, 0.95, 0.95]', ...
   'FontWeight', 'bold');
note = Simulink.Annotation('xmat_test', 'Output_values');
set(note, 'FixedHeight', 'on', 'FixedWidth', 'on', 'Position', [665, 85, 835, 535], ...
   'ForegroundColor', '[0.45, 0.45, 0.45]', 'BackgroundColor', '[0.95, 0.95, 0.95]', ...
   'FontWeight', 'bold');

save_system('xmat_test', 'xmat_test.slx')