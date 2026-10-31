%% STEP 27.8 - Create Five-Cell Simulink Thermal Model
% Battery Thermal Safety Assessment: LFP vs NMC811
%
% Programmatically builds a Simulink model using a
% five-state linear thermal network.
%
% Initial implementation: LFP, constant 100 W heater.
% No thermal runaway reaction included.

clc;

%% 1. Load parameters

step01_thermal_model_parameters;

battery = LFP;

n = model.nCells;

C = battery.thermal_capacity_J_K;
R = model.R_intercell_K_W;
h = model.h_W_m2K;
Tamb = model.T_ambient_C;

A_surface = 2 * ( ...
    battery.length_m * battery.width_m + ...
    battery.length_m * battery.height_m + ...
    battery.width_m * battery.height_m);

%% 2. Build state-space matrices

% State vector: temperature rise above ambient [K]
% x = T - Tamb
%
% C*dx/dt = -L*x - h*A*x + B*Q

A_matrix = zeros(n,n);

for i = 1:n

    neighbor_count = 0;

    if i > 1
        A_matrix(i,i-1) = 1 / (R*C);
        neighbor_count = neighbor_count + 1;
    end

    if i < n
        A_matrix(i,i+1) = 1 / (R*C);
        neighbor_count = neighbor_count + 1;
    end

    A_matrix(i,i) = ...
        -(neighbor_count/R + h*A_surface) / C;

end

B_matrix = zeros(n,1);
B_matrix(1) = 1/C;

C_matrix = eye(n);
D_matrix = zeros(n,1);

%% 3. Create Simulink model

model_name = 'five_cell_thermal_network';

model_file = fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))), ...
    'models', ...
    [model_name '.slx']);

if bdIsLoaded(model_name)
    close_system(model_name, 0);
end

new_system(model_name);
open_system(model_name);

%% 4. Add blocks

add_block( ...
    'simulink/Sources/Constant', ...
    [model_name '/Heater'], ...
    'Value', '100', ...
    'Position', [80 110 130 140]);

add_block( ...
    'simulink/Continuous/State-Space', ...
    [model_name '/Thermal Network'], ...
    'A', mat2str(A_matrix), ...
    'B', mat2str(B_matrix), ...
    'C', mat2str(C_matrix), ...
    'D', mat2str(D_matrix), ...
    'X0', 'zeros(5,1)', ...
    'Position', [210 90 370 160]);

add_block( ...
    'simulink/Math Operations/Bias', ...
    [model_name '/Add Ambient Temperature'], ...
    'Bias', num2str(Tamb), ...
    'Position', [440 105 560 145]);

add_block( ...
    'simulink/Sinks/Scope', ...
    [model_name '/Temperature Scope'], ...
    'Position', [650 105 700 145]);

add_block( ...
    'simulink/Sinks/To Workspace', ...
    [model_name '/Temperature Data'], ...
    'VariableName', 'sim_temperature', ...
    'SaveFormat', 'Timeseries', ...
    'Position', [650 205 760 235]);

%% 5. Connect blocks

add_line(model_name, ...
    'Heater/1', 'Thermal Network/1');

add_line(model_name, ...
    'Thermal Network/1', ...
    'Add Ambient Temperature/1');

add_line(model_name, ...
    'Add Ambient Temperature/1', ...
    'Temperature Scope/1');

add_line(model_name, ...
    'Add Ambient Temperature/1', ...
    'Temperature Data/1');

%% 6. Configure simulation

set_param(model_name, ...
    'StopTime', '2000', ...
    'Solver', 'ode45');

%% 7. Save model

save_system(model_name, model_file);

fprintf('\nSIMULINK MODEL CREATED\n');
fprintf('Model: %s\n', model_name);
fprintf('Saved to: %s\n', model_file);
fprintf('\nSTEP 27.8 COMPLETED\n');
