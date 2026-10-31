%% STEP 27.19 - Simplified Thermal Runaway Model
% Battery Thermal Safety Assessment: LFP vs NMC811
%
% Five-cell lumped thermal network with triggered heat release.
%
% IMPORTANT:
% This is an illustrative phenomenological model.
% TR parameters are NOT experimentally calibrated.
% The model does not represent detailed reaction kinetics.

clc;

%% 1. Load battery parameters

step01_thermal_model_parameters;

%% 2. Select battery chemistry

battery = LFP;

% To simulate NMC811 instead, change to:
% battery = NMC;

%% 3. Simulation settings

t_end = 3000;
Q_heater_W = 100;

T_ambient = model.T_ambient_C;
T_initial = model.T_initial_C;

n = model.nCells;

%% 4. Thermal properties

C = battery.mass_kg * battery.cp_J_kgK;

R = model.R_intercell_K_W;
h = model.h_W_m2K;

A = 2 * ( ...
    battery.length_m * battery.width_m + ...
    battery.length_m * battery.height_m + ...
    battery.width_m * battery.height_m);

%% 5. Illustrative Thermal Runaway parameters

TR.T_trigger_C = 70;

TR.Q_max_W = 500;

TR.duration_s = 300;

% Total heat released per triggered cell:
TR.E_release_J = TR.Q_max_W * TR.duration_s;

fprintf('\n========================================\n');
fprintf('SIMPLIFIED THERMAL RUNAWAY MODEL\n');
fprintf('========================================\n');

fprintf('Chemistry: %s\n', battery.name);
fprintf('Trigger temperature: %.1f C\n', TR.T_trigger_C);
fprintf('TR heat release power: %.1f W\n', TR.Q_max_W);
fprintf('TR duration: %.1f s\n', TR.duration_s);
fprintf('TR energy per cell: %.1f kJ\n', ...
    TR.E_release_J / 1000);

%% 6. Initial conditions

% Temperature states
T0 = T_initial * ones(n, 1);

% Released-energy states
E0 = zeros(n, 1);

% Combined state vector:
% y = [T1 ... T5 E1 ... E5]
y0 = [T0; E0];

%% 7. Solve ODE

% Limit maximum step to improve trigger resolution
options = odeset('MaxStep', 0.1, ...
                 'RelTol', 1e-6, ...
                 'AbsTol', 1e-8);

[t, y] = ode45( ...
    @(t, y) thermal_runaway_ode( ...
        t, y, C, R, h, A, ...
        T_ambient, Q_heater_W, TR, n), ...
    [0 t_end], y0, options);

%% 8. Extract temperatures

T = y(:, 1:n);

%% 9. Plot temperatures

figure('Color', 'w', ...
    'Name', 'Simplified Thermal Runaway');

plot(t, T, 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Cell temperature [°C]');

title(['Simplified Thermal Runaway - ' battery.name]);

legend('Cell 1', 'Cell 2', 'Cell 3', ...
       'Cell 4', 'Cell 5', ...
       'Location', 'best');

grid on;
box on;
%% 9.1 Verify Thermal Runaway energy release

E = y(:, n+1:2*n);

fprintf('\n========================================\n');
fprintf('THERMAL RUNAWAY ENERGY VERIFICATION\n');
fprintf('========================================\n');

fprintf('Maximum allowed energy: %.2f kJ per cell\n', ...
    TR.E_release_J / 1000);

for i = 1:n

    E_released_kJ = E(end, i) / 1000;

    fprintf('Cell %d released energy: %.4f kJ\n', ...
        i, E_released_kJ);

end
%% 10. Display results

fprintf('\nFinal cell temperatures:\n');

for i = 1:n
    fprintf('Cell %d: %.2f C\n', i, T(end, i));
end

fprintf('\nMaximum cell temperatures:\n');

for i = 1:n
    fprintf('Cell %d: %.2f C\n', i, max(T(:, i)));
end

fprintf('\nSTEP 27.19 COMPLETED\n');

%% Local ODE function

function dydt = thermal_runaway_ode( ...
    ~, y, C, R, h, A, Tamb, Qheater, TR, n)

    T = y(1:n);
    E = y(n+1:2*n);

    dTdt = zeros(n, 1);
    dEdt = zeros(n, 1);

    for i = 1:n

        %% External heating

        Q_external = 0;

        if i == 1
            Q_external = Qheater;
        end

        %% Thermal Runaway heat release

        Q_TR = 0;

        if T(i) >= TR.T_trigger_C && ...
                E(i) < TR.E_release_J

            Q_TR = TR.Q_max_W;

        end

        %% Neighbor heat transfer

        Q_neighbors = 0;

        if i > 1
            Q_neighbors = Q_neighbors + ...
                (T(i-1) - T(i)) / R;
        end

        if i < n
            Q_neighbors = Q_neighbors + ...
                (T(i+1) - T(i)) / R;
        end

        %% Ambient heat losses

        Q_loss = h * A * (T(i) - Tamb);

        %% Energy balance

        dTdt(i) = ( ...
            Q_external + Q_TR + ...
            Q_neighbors - Q_loss) / C;

        %% Accumulated TR energy

        dEdt(i) = Q_TR;

    end

    dydt = [dTdt; dEdt];

end
