%% STEP 27.5 - Five-Cell Thermal Network
% Battery Thermal Safety Assessment: LFP vs NMC811
%
% Simplified lumped thermal model.
% No thermal runaway reaction is included at this stage.

clc;

%% 1. Load model parameters

step01_thermal_model_parameters;

%% 2. Select chemistry

battery = NMC;

%% 3. Simulation settings

tspan = [0 2000];

T0 = model.T_initial_C * ones(model.nCells, 1);

% Illustrative external heating of Cell 1
Q_heater_W = 100;

%% 4. Calculate cell surface area

A = 2 * ( ...
    battery.length_m * battery.width_m + ...
    battery.length_m * battery.height_m + ...
    battery.width_m * battery.height_m);

%% 5. Define thermal model

C = battery.thermal_capacity_J_K;
R = model.R_intercell_K_W;
h = model.h_W_m2K;
Tamb = model.T_ambient_C;

thermalODE = @(t, T) thermal_network_ode( ...
    t, T, C, R, h, A, Tamb, Q_heater_W);

%% 6. Solve differential equations

[t, T] = ode45(thermalODE, tspan, T0);

%% 7. Plot results

figure;

plot(t, T, 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Temperature [°C]');
title(['Five-Cell Thermal Network - ' battery.name]);

legend('Cell 1', 'Cell 2', 'Cell 3', ...
       'Cell 4', 'Cell 5', 'Location', 'best');

grid on;

%% 8. Display final temperatures

fprintf('\nFIVE-CELL THERMAL SIMULATION\n');

for i = 1:model.nCells
    fprintf('Cell %d: %.2f °C\n', i, T(end,i));
end
%% 9. Energy balance verification

% Total energy supplied by the heater
E_heater_J = Q_heater_W * (t(end) - t(1));

% Increase in stored thermal energy
E_stored_J = C * sum(T(end,:) - T(1,:));

% Instantaneous heat loss to ambient for all cells
Q_ambient_W = h * A * sum(T - Tamb, 2);

% Numerically integrated heat loss
E_ambient_J = trapz(t, Q_ambient_W);

% Energy balance residual
E_residual_J = E_heater_J - E_stored_J - E_ambient_J;

relative_error_percent = ...
    abs(E_residual_J) / E_heater_J * 100;

fprintf('\nENERGY BALANCE VERIFICATION\n');
fprintf('Heater energy: %.2f J\n', E_heater_J);
fprintf('Stored energy: %.2f J\n', E_stored_J);
fprintf('Ambient losses: %.2f J\n', E_ambient_J);
fprintf('Residual: %.2f J\n', E_residual_J);
fprintf('Relative error: %.6f %%\n', relative_error_percent);
fprintf('\nSTEP 27.5 COMPLETED\n');

%% Local function

function dTdt = thermal_network_ode( ...
    ~, T, C, R, h, A, Tamb, Qheater)

    n = length(T);

    dTdt = zeros(n,1);

    for i = 1:n

        % External heating
        Q_in = 0;

        if i == 1
            Q_in = Qheater;
        end

        % Heat transfer between adjacent cells
        Q_neighbors = 0;

        if i > 1
            Q_neighbors = Q_neighbors + ...
                (T(i-1) - T(i)) / R;
        end

        if i < n
            Q_neighbors = Q_neighbors + ...
                (T(i+1) - T(i)) / R;
        end

        % Heat loss to ambient
        Q_loss = h * A * (T(i) - Tamb);

        % Energy balance
        dTdt(i) = (Q_in + Q_neighbors - Q_loss) / C;

    end
end