%% STEP 27 FINAL - Event-driven thermal runaway propagation demonstration
% Five-cell lumped thermal model, illustrative LFP geometry.
% NOT calibrated to experimental TR kinetics or onset temperatures.
% Two scenarios differ in intercell resistance only.
clc;
step01_thermal_model_parameters;
battery = LFP;
n = model.nCells;
Tamb = model.T_ambient_C;
C = battery.mass_kg * battery.cp_J_kgK;
A = 2 * (battery.length_m*battery.width_m + ...
         battery.length_m*battery.height_m + ...
         battery.width_m*battery.height_m);
h = model.h_W_m2K;
TR.trigger_C = 65;               % Illustrative, NOT experimental
TR.power_W = 1000;                % Illustrative, NOT experimental
TR.energy_J = 500e3;             % Illustrative per-cell budget
heater_W = 100;                  % Applied continuously to cell 1
tEnd = 6000;
% Lower R allows more energy transfer. Both are illustrative.
scenarios = struct('name', {'Isolated', 'Coupled'}, ...
                   'R', {2.0, 0.5});
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(fileparts(scriptDir), 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
summary = table();
fig = figure('Color', 'w', 'Name', 'TR propagation scenarios');
colors = lines(n);
for s = 1:numel(scenarios)
    R = scenarios(s).R;
    tAll = 0;
    y = [model.T_initial_C * ones(n,1); zeros(n,1)];
    yAll = y.';
    % mode: 0 inactive, 1 releasing, 2 exhausted
    mode = zeros(n,1);
    eventTime = nan(n,1);
    % The initial state might be above trigger (not the case here).
    for k = 1:500
        t0 = tAll(end);
        if t0 >= tEnd - 1e-9, break; end
        % Detect only relevant event surfaces; preserve event->cell mapping.
        activeIdx = find(mode == 1);
        inactiveIdx = find(mode == 0);
        eventCells = [inactiveIdx; activeIdx];
        eventKinds = [ones(numel(inactiveIdx),1); 2*ones(numel(activeIdx),1)];
        opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, ...
            'MaxStep', 1, 'Events', @(t,z) events(t,z,eventCells,eventKinds,TR,n));
        [tt, yy, te, ye, ie] = ode45(@(t,z) rhs(z,mode,C,R,h,A,Tamb,heater_W,TR,n), ...
                                    [t0 tEnd], y, opts);
        if numel(tt)>1
            tAll = [tAll; tt(2:end)]; %#ok<AGROW>
            yAll = [yAll; yy(2:end,:)]; %#ok<AGROW>
        end
        y = yy(end,:).';
        if isempty(ie), break; end
        % Process all events at the same time; update state at exact surface.
        hit = unique(ie(:));
        for q = 1:numel(hit)
            cellIdx = eventCells(hit(q));
            if eventKinds(hit(q)) == 1 && mode(cellIdx) == 0
                mode(cellIdx) = 1;
                eventTime(cellIdx) = tt(end);
            elseif eventKinds(hit(q)) == 2 && mode(cellIdx) == 1
                mode(cellIdx) = 2;
                y(n+cellIdx) = TR.energy_J;
            end
        end
        % Catch any simultaneously crossed thresholds after a discontinuity.
        newlyHot = find(mode == 0 & y(1:n) >= TR.trigger_C - 1e-7);
        for q = 1:numel(newlyHot)
            i = newlyHot(q);
            mode(i) = 1;
            eventTime(i) = tt(end);
        end
        exhausted = find(mode == 1 & y(n+1:end) >= TR.energy_J - 1e-5);
        for q = 1:numel(exhausted)
            i = exhausted(q);
            mode(i) = 2;
            y(n+i) = TR.energy_J;
        end
        % Update the last saved row to reflect event-projected energy state.
        yAll(end,:) = y.';
    end
    if tAll(end) < tEnd - 1e-6
        error('Scenario %s did not finish; inspect event iteration.', scenarios(s).name);
    end
    for i = 1:n
        summary = [summary; table(string(scenarios(s).name), R, i, ...
            eventTime(i), max(yAll(:,i)), yAll(end,n+i)/1000, ...
            'VariableNames', {'Scenario','R_K_per_W','Cell', ...
            'TriggerTime_s','MaxTemperature_C','ReleasedEnergy_kJ'})]; %#ok<AGROW>
    end
    subplot(1,2,s);
    hold on;
    for i = 1:n
        plot(tAll, yAll(:,i), 'LineWidth', 1.4, 'Color', colors(i,:));
    end
    yline(TR.trigger_C, 'k--', 'Trigger threshold');
    xlabel('Time [s]'); ylabel('Temperature [°C]');
    title(sprintf('%s: R = %.3g K/W', scenarios(s).name, R));
    grid on; box on;
    if s==2, legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5', 'Location','best'); end
    fprintf('\n%s (R=%.3g K/W): %d/%d cells triggered.\n', ...
        scenarios(s).name, R, sum(~isnan(eventTime)), n);
    fprintf('Released energies [kJ]: ');
    fprintf('%.4f ', yAll(end,n+1:end)/1000);
    fprintf('\n');
end
writetable(summary, fullfile(resultsDir, 'thermal_runaway_scenarios.csv'));
exportgraphics(fig, fullfile(resultsDir, 'thermal_runaway_scenarios.png'), 'Resolution', 300);
disp(summary);
fprintf('\nSTEP 27 COMPLETE: results exported to %s\n', resultsDir);

function dz = rhs(z,mode,C,R,h,A,Tamb,Qheater,TR,n)
    T = z(1:n);
    Q = zeros(n,1); Q(1) = Qheater;
    Q = Q + TR.power_W * (mode == 1);
    Q = Q - h*A*(T-Tamb);
    for i=1:n-1
        flux = (T(i)-T(i+1))/R;
        Q(i) = Q(i)-flux;
        Q(i+1) = Q(i+1)+flux;
    end
    dE = TR.power_W * (mode == 1);
    dz = [Q/C; dE];
end

function [value,isterminal,direction] = events(~,z,cells,kinds,TR,n)
    value = zeros(numel(cells),1);
    direction = zeros(numel(cells),1);
    for j=1:numel(cells)
        i = cells(j);
        if kinds(j)==1
            value(j) = z(i)-TR.trigger_C;
            direction(j) = +1;
        else
            value(j) = z(n+i)-TR.energy_J;
            direction(j) = +1;
        end
    end
    isterminal = ones(size(value));
end
