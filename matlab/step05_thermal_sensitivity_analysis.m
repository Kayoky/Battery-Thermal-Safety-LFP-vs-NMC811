%% STEP 27.17 - Thermal Sensitivity Analysis
% Battery Thermal Safety Assessment: LFP vs NMC811
%
% Five-cell lumped thermal network
% One-factor-at-a-time sensitivity analysis
%
% IMPORTANT:
% Thermal parameters are illustrative assumptions.
% This is NOT a validated thermal runaway model.

clc;

%% 1. Load parameters

step01_thermal_model_parameters;

%% 2. Simulation settings

tspan = [0 2000];
Q_heater_W = 100;

batteries = {LFP, NMC};

% Parameter ranges
cp_values = [700 1000 1300];       % J/(kg*K)
R_values = [0.5 2.0 5.0];          % K/W
h_values = [5 10 20];              % W/(m^2*K)

% Baseline parameters
cp_base = 1000;
R_base = 2.0;
h_base = 10;

%% 3. Output directory

script_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(script_dir);
results_dir = fullfile(matlab_dir, 'results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% 4. Run sensitivity analysis

% Store results in a table
sensitivity_results = table();

for b = 1:numel(batteries)

    battery = batteries{b};

    % Geometry
    A = 2 * ( ...
        battery.length_m * battery.width_m + ...
        battery.length_m * battery.height_m + ...
        battery.width_m * battery.height_m);

    % Define parameter studies
    studies = {
        'Specific heat capacity', cp_values;
        'Intercell resistance', R_values;
        'Heat transfer coefficient', h_values
    };

    for s = 1:size(studies, 1)

        parameter_name = studies{s, 1};
        parameter_values = studies{s, 2};

        for p = 1:numel(parameter_values)

            % Reset to baseline for every run
            cp = cp_base;
            R = R_base;
            h = h_base;

            value = parameter_values(p);

            % Change one parameter at a time
            switch parameter_name

                case 'Specific heat capacity'
                    cp = value;

                case 'Intercell resistance'
                    R = value;

                case 'Heat transfer coefficient'
                    h = value;

            end

            % Thermal capacity
            C = battery.mass_kg * cp;

            % Initial temperatures
            T0 = model.T_initial_C * ...
                ones(model.nCells, 1);

            % Thermal network ODE
            ode_function = @(t, T) thermal_ode( ...
                T, C, R, h, A, ...
                model.T_ambient_C, Q_heater_W);

            % Solve
            [t, T] = ode45(ode_function, tspan, T0);

            % Extract final temperatures
            T1_final = T(end, 1);
            T2_final = T(end, 2);
            T5_final = T(end, 5);

            % Store result
            new_row = table( ...
                string(battery.name), ...
                string(parameter_name), ...
                value, ...
                T1_final, ...
                T2_final, ...
                T5_final, ...
                'VariableNames', { ...
                'Chemistry', ...
                'Parameter', ...
                'Value', ...
                'Cell1_C', ...
                'Cell2_C', ...
                'Cell5_C'});

            sensitivity_results = [
                sensitivity_results;
                new_row
            ];

        end
    end
end

%% 5. Save CSV

csv_file = fullfile( ...
    results_dir, 'thermal_sensitivity_results.csv');

writetable(sensitivity_results, csv_file);

%% 6. Display results

fprintf('\n========================================\n');
fprintf('THERMAL SENSITIVITY ANALYSIS\n');
fprintf('========================================\n');

fprintf('Total simulations: %d\n', ...
    height(sensitivity_results));

disp(sensitivity_results);

fprintf('\nResults saved to:\n%s\n', csv_file);

fprintf('\nSTEP 27.17 COMPLETED\n');

%% Local thermal model

function dTdt = thermal_ode( ...
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