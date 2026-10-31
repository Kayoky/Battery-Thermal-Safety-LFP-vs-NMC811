%% STEP 27.3 - Thermal Model Parameters
% Battery Thermal Safety Assessment: LFP vs NMC811
% Five-cell lumped thermal network
%
% Experimental cell specifications:
% Schoberl et al., eTransportation 19 (2024), 100305
%
% IMPORTANT:
% Thermal parameters below are preliminary modeling assumptions.
% They are NOT experimentally identified material properties.

clc;

%% 1. General model settings

model.nCells = 5;
model.T_ambient_C = 25;
model.T_initial_C = 25;

%% 2. Experimental cell specifications

% LFP
LFP.name = 'LFP';
LFP.mass_kg = 3.10;
LFP.length_m = 0.280;
LFP.width_m = 0.063;
LFP.height_m = 0.082;
LFP.nominal_voltage_V = 3.20;
LFP.capacity_Ah = 163.3;

% NMC811
NMC.name = 'NMC811';
NMC.mass_kg = 1.42;
NMC.length_m = 0.173;
NMC.width_m = 0.042;
NMC.height_m = 0.085;
NMC.nominal_voltage_V = 3.65;
NMC.capacity_Ah = 96.3;

%% 3. Preliminary thermal assumptions

% Illustrative values for initial model development.
% These must be varied in sensitivity studies.

LFP.cp_J_kgK = 1000;
NMC.cp_J_kgK = 1000;

model.h_W_m2K = 10;

% Effective thermal resistance between adjacent cells
model.R_intercell_K_W = 2.0;

%% 4. Derived parameters

LFP.thermal_capacity_J_K = LFP.mass_kg * LFP.cp_J_kgK;
NMC.thermal_capacity_J_K = NMC.mass_kg * NMC.cp_J_kgK;

LFP.volume_m3 = LFP.length_m * LFP.width_m * LFP.height_m;
NMC.volume_m3 = NMC.length_m * NMC.width_m * NMC.height_m;

%% 5. Display parameters

fprintf('\n========================================\n');
fprintf('FIVE-CELL THERMAL MODEL PARAMETERS\n');
fprintf('========================================\n');

fprintf('\nLFP:\n');
fprintf('Mass: %.2f kg\n', LFP.mass_kg);
fprintf('Thermal capacity: %.1f J/K\n', ...
    LFP.thermal_capacity_J_K);
fprintf('Volume: %.6f m^3\n', LFP.volume_m3);

fprintf('\nNMC811:\n');
fprintf('Mass: %.2f kg\n', NMC.mass_kg);
fprintf('Thermal capacity: %.1f J/K\n', ...
    NMC.thermal_capacity_J_K);
fprintf('Volume: %.6f m^3\n', NMC.volume_m3);

fprintf('\nModel assumptions:\n');
fprintf('Ambient temperature: %.1f C\n', ...
    model.T_ambient_C);
fprintf('Heat transfer coefficient: %.1f W/(m^2 K)\n', ...
    model.h_W_m2K);
fprintf('Inter-cell thermal resistance: %.2f K/W\n', ...
    model.R_intercell_K_W);

fprintf('\nSTEP 27.3 COMPLETED\n');
