%% STEP 27.18 - Thermal Sensitivity Visualization
% Battery Thermal Safety Assessment: LFP vs NMC811
%
% Visualizes the one-factor-at-a-time sensitivity analysis.
% All results are based on illustrative thermal parameters.
% This is NOT a validated thermal runaway model.

clc;

%% 1. Define paths

script_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(script_dir);
results_dir = fullfile(matlab_dir, 'results');

csv_file = fullfile(results_dir, ...
    'thermal_sensitivity_results.csv');

%% 2. Load sensitivity results

if ~isfile(csv_file)
    error('Sensitivity results not found: %s', csv_file);
end

data = readtable(csv_file, 'TextType', 'string');

fprintf('\n========================================\n');
fprintf('THERMAL SENSITIVITY VISUALIZATION\n');
fprintf('========================================\n');

fprintf('Loaded %d simulation results.\n', height(data));

%% 3. Define parameter studies

parameters = [
    "Specific heat capacity"
    "Intercell resistance"
    "Heat transfer coefficient"
];

x_labels = [
    "Specific heat capacity [J/(kg K)]"
    "Intercell thermal resistance [K/W]"
    "Heat transfer coefficient [W/(m^2 K)]"
];

file_names = [
    "sensitivity_specific_heat_capacity.png"
    "sensitivity_intercell_resistance.png"
    "sensitivity_heat_transfer_coefficient.png"
];

%% 4. Generate figures

for p = 1:numel(parameters)

    parameter_name = parameters(p);

    % Select results for this parameter
    selected = data(data.Parameter == parameter_name, :);

    % Separate chemistries
    LFP_data = selected(selected.Chemistry == "LFP", :);
    NMC_data = selected(selected.Chemistry == "NMC811", :);

    % Sort by parameter value
    LFP_data = sortrows(LFP_data, 'Value');
    NMC_data = sortrows(NMC_data, 'Value');

    % Create figure
    fig = figure('Color', 'w', ...
        'Name', char(parameter_name));

    plot(LFP_data.Value, LFP_data.Cell1_C, ...
        '-o', 'LineWidth', 1.8, 'MarkerSize', 7);

    hold on;

    plot(NMC_data.Value, NMC_data.Cell1_C, ...
        '-s', 'LineWidth', 1.8, 'MarkerSize', 7);

    xlabel(x_labels(p));
    ylabel('Cell 1 final temperature [°C]');

    title("Thermal Sensitivity: " + parameter_name);

    legend('LFP', 'NMC811', 'Location', 'best');

    grid on;
    box on;

    % Display actual tested parameter values
    xticks(unique(selected.Value));

    % Export figure
    output_file = fullfile(results_dir, file_names(p));

    exportgraphics(fig, output_file, 'Resolution', 300);

    fprintf('Saved: %s\n', file_names(p));

end

fprintf('\nSTEP 27.18 COMPLETED\n');