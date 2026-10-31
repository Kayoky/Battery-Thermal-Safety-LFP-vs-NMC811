%% STEP 27.15 - LFP vs NMC811 Thermal Comparison
% Battery Thermal Safety Assessment
%
% Five-cell lumped thermal network
%
% This simulation compares the two experimental cell formats
% under identical illustrative thermal boundary conditions.
%
% IMPORTANT:
% This is NOT a validated thermal runaway simulation.
% Differences reflect cell mass and geometry as well as
% the assumptions used in the thermal model.

clc;

%% 1. Load parameters

step01_thermal_model_parameters;

%% 2. Simulation settings

tspan = [0 2000];

Q_heater_W = 100;

batteries = {LFP, NMC};

results = struct();

%% 3. Run simulations

for k = 1:numel(batteries)

    battery = batteries{k};

    n = model.nCells;

    C = battery.thermal_capacity_J_K;
    R = model.R_intercell_K_W;
    h = model.h_W_m2K;
    Tamb = model.T_ambient_C;

    % Total external surface area of one cell
    A = 2 * ( ...
        battery.length_m * battery.width_m + ...
        battery.length_m * battery.height_m + ...
        battery.width_m * battery.height_m);

    % Initial temperatures
    T0 = model.T_initial_C * ones(n, 1);

    % Thermal network equations
    thermalODE = @(t, T) local_thermal_ode( ...
        T, C, R, h, A, Tamb, Q_heater_W);

    % Numerical integration
    [t, T] = ode45(thermalODE, tspan, T0);

    % Store results
    results.(battery.name).time = t;
    results.(battery.name).temperature = T;
    results.(battery.name).final_temperature = T(end, :);

end

%% 4. Display comparison

fprintf('\n========================================\n');
fprintf('LFP VS NMC811 THERMAL COMPARISON\n');
fprintf('========================================\n');

fprintf('\nHeater power: %.1f W\n', Q_heater_W);
fprintf('Simulation time: %.1f s\n', tspan(2));

fprintf('\nCell       LFP [C]       NMC811 [C]\n');
fprintf('--------------------------------------\n');

for i = 1:model.nCells

    fprintf('%d          %8.2f       %8.2f\n', ...
        i, ...
        results.LFP.final_temperature(i), ...
        results.NMC811.final_temperature(i));

end

%% 5. Plot Cell 1 comparison

figure('Name', 'LFP vs NMC811 - Cell 1');

plot( ...
    results.LFP.time, ...
    results.LFP.temperature(:,1), ...
    'LineWidth', 1.8);

hold on;

plot( ...
    results.NMC811.time, ...
    results.NMC811.temperature(:,1), ...
    'LineWidth', 1.8);

xlabel('Time [s]');
ylabel('Temperature [°C]');

title('Cell 1 Thermal Response: LFP vs NMC811');

legend('LFP', 'NMC811', 'Location', 'best');

grid on;

%% 6. Plot final temperatures

figure('Name', 'Final Cell Temperatures');

cell_numbers = 1:model.nCells;

plot( ...
    cell_numbers, ...
    results.LFP.final_temperature, ...
    '-o', 'LineWidth', 1.8);

hold on;

plot( ...
    cell_numbers, ...
    results.NMC811.final_temperature, ...
    '-s', 'LineWidth', 1.8);

xlabel('Cell Number');
ylabel('Final Temperature [°C]');

title('Final Cell Temperatures After 2000 s');

legend('LFP', 'NMC811', 'Location', 'best');

xticks(cell_numbers);

grid on;
%% 6.1 Save comparison figures

script_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(script_dir);
results_dir = fullfile(matlab_dir, 'results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Find figures
fig_cell1 = findobj('Type', 'figure', ...
    'Name', 'LFP vs NMC811 - Cell 1');

fig_final = findobj('Type', 'figure', ...
    'Name', 'Final Cell Temperatures');

% Export figures
exportgraphics(fig_cell1, fullfile(results_dir, ...
    'thermal_cell1_comparison.png'), 'Resolution', 300);

exportgraphics(fig_final, fullfile(results_dir, ...
    'thermal_final_temperature_comparison.png'), 'Resolution', 300);

fprintf('\nComparison figures exported successfully.\n');
%% 7. Save numerical results

script_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(script_dir);
results_dir = fullfile(matlab_dir, 'results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Use common time points for the exported temperature curves
t_common = (0:10:2000)';

LFP_interp = interp1( ...
    results.LFP.time, ...
    results.LFP.temperature, ...
    t_common, 'linear');

NMC_interp = interp1( ...
    results.NMC811.time, ...
    results.NMC811.temperature, ...
    t_common, 'linear');

% Export time histories
column_names = [ ...
    "Time_s", ...
    "LFP_Cell1_C", "LFP_Cell2_C", "LFP_Cell3_C", ...
    "LFP_Cell4_C", "LFP_Cell5_C", ...
    "NMC811_Cell1_C", "NMC811_Cell2_C", ...
    "NMC811_Cell3_C", "NMC811_Cell4_C", ...
    "NMC811_Cell5_C"];

comparison_table = array2table( ...
    [t_common, LFP_interp, NMC_interp], ...
    'VariableNames', cellstr(column_names));

writetable( ...
    comparison_table, ...
    fullfile(results_dir, 'thermal_comparison.csv'));

% Export final temperatures
final_table = table( ...
    cell_numbers', ...
    results.LFP.final_temperature', ...
    results.NMC811.final_temperature', ...
    'VariableNames', ...
    {'Cell', 'LFP_C', 'NMC811_C'});

writetable( ...
    final_table, ...
    fullfile(results_dir, 'thermal_final_temperatures.csv'));

fprintf('\nResults saved to:\n%s\n', results_dir);

fprintf('\nSTEP 27.15 COMPLETED\n');

%% Local thermal model function

function dTdt = local_thermal_ode( ...
    T, C, R, h, A, Tamb, Qheater)

    n = numel(T);
    dTdt = zeros(n, 1);

    for i = 1:n

        Q_in = 0;

        if i == 1
            Q_in = Qheater;
        end

        Q_neighbors = 0;

        if i > 1
            Q_neighbors = Q_neighbors + ...
                (T(i-1) - T(i)) / R;
        end

        if i < n
            Q_neighbors = Q_neighbors + ...
                (T(i+1) - T(i)) / R;
        end

        Q_loss = h * A * (T(i) - Tamb);

        dTdt(i) = ...
            (Q_in + Q_neighbors - Q_loss) / C;

    end

end