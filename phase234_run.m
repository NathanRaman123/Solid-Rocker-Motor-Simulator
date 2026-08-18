%% phase234_run.m — Phases 2, 3, 4 Master Demo Script
%
%  Runs the full simulation pipeline:
%    Phase 2: Propellant thermochemistry + nozzle model
%    Phase 3: Internal ballistics ODE (ode45)
%    Phase 4: Trade studies, validation, results comparison
%
%  Prerequisites: All .m files from Phase 1 and this phase in MATLAB path.
%  Toolboxes: Optimization Toolbox (for fsolve). If unavailable, see note
%             at bottom of this file for a fsolve substitute.
%
%  Run time: ~10–30 seconds depending on sweep size.
% ─────────────────────────────────────────────────────────────────────────

clear; clc; close all;
fprintf('=== Solid Rocket Motor Simulator — Phases 2/3/4 ===\n\n');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 2 — Propellant thermochemistry model
%  ════════════════════════════════════════════════════════════════════════

fprintf('--- PHASE 2: Propellant Model ---\n\n');

% Load two propellants for comparison
prop_apcp = propellant_db('APCP_Aerotech');
prop_knsb = propellant_db('KNSB');

% Burn rate curves for both
fprintf('\n[Burn rate at 3.5 MPa]\n');
r_apcp = burn_rate(prop_apcp, 3.5e6);
r_knsb = burn_rate(prop_knsb, 3.5e6);

% Plot burn rate curves with temperature sensitivity bands
figure('Name', 'Phase 2: Propellant Comparison', 'Color', [1 1 1], ...
       'Position', [50 50 1000 420]);

subplot(1,2,1);
[~, rc_apcp] = burn_rate(prop_apcp, []);
hold on; grid on; box on;
fill([rc_apcp.P, fliplr(rc_apcp.P)]/1e6, ...
     [rc_apcp.r_cold, fliplr(rc_apcp.r_hot)]*1e3, ...
     [0.8 0.88 1.0], 'EdgeColor','none','FaceAlpha',0.5,'DisplayName','APCP ±30K band');
plot(rc_apcp.P/1e6, rc_apcp.r*1e3, 'b-', 'LineWidth', 2.5, 'DisplayName', 'APCP Aerotech');

[~, rc_knsb] = burn_rate(prop_knsb, []);
fill([rc_knsb.P, fliplr(rc_knsb.P)]/1e6, ...
     [rc_knsb.r_cold, fliplr(rc_knsb.r_hot)]*1e3, ...
     [1.0 0.88 0.80], 'EdgeColor','none','FaceAlpha',0.5,'DisplayName','KNSB ±30K band');
plot(rc_knsb.P/1e6, rc_knsb.r*1e3, 'r-', 'LineWidth', 2.5, 'DisplayName', 'KNSB');

xlabel('Chamber pressure [MPa]'); ylabel('Burn rate [mm/s]');
title("Saint-Robert's law: r = a·P^n", 'FontWeight','normal');
legend('Location','northwest'); xlim([0 12]);

% Nozzle model for both propellants
At      = pi * (0.015/2)^2;   % 15mm throat diameter
epsilon = 8;                   % expansion ratio

nozzle_apcp = nozzle_model(prop_apcp, At, epsilon, 101325);
nozzle_knsb = nozzle_model(prop_knsb, At, epsilon, 101325);

subplot(1,2,2);
Pc_range = linspace(0.5e6, 8e6, 200);
Cf_apcp  = arrayfun(@(P) nozzle_apcp.Cf(P), Pc_range);
Cf_knsb  = arrayfun(@(P) nozzle_knsb.Cf(P), Pc_range);
hold on; grid on; box on;
plot(Pc_range/1e6, Cf_apcp, 'b-', 'LineWidth', 2.5, 'DisplayName', 'APCP nozzle');
plot(Pc_range/1e6, Cf_knsb, 'r-', 'LineWidth', 2.5, 'DisplayName', 'KNSB nozzle');
xlabel('Chamber pressure [MPa]'); ylabel('Thrust coefficient C_f');
title(sprintf('Thrust coefficient (ε=%.0f, Dt=%.0fmm)', epsilon, nozzle_apcp.Dt*1e3), ...
      'FontWeight','normal');
legend('Location','southeast');
sgtitle('Phase 2 — Propellant & Nozzle Characterization', ...
        'FontSize',12,'FontWeight','normal');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 3 — Internal ballistics ODE
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n--- PHASE 3: Internal Ballistics ODE ---\n\n');

% Define grain (BATES, 4-segment)
g_bates = BATES(0.075, 0.025, 0.150, 4, 500);

% Simulate: APCP Aerotech
res_apcp = ballistics_ode(g_bates, prop_apcp, nozzle_apcp, 'plot', true);

% Simulate: KNSB (same grain, different propellant)
res_knsb = ballistics_ode(g_bates, prop_knsb, nozzle_knsb);

% Side-by-side thrust comparison
figure('Name','Phase 3: APCP vs KNSB — Same Grain', 'Color',[1 1 1], ...
       'Position',[100 150 900 380]);
subplot(1,2,1);
hold on; grid on; box on;
plot(res_apcp.t, res_apcp.F, 'b-', 'LineWidth', 2, ...
     'DisplayName', sprintf('APCP | It=%.0fNs | Class %s', res_apcp.It, res_apcp.motor_class));
plot(res_knsb.t, res_knsb.F, 'r-', 'LineWidth', 2, ...
     'DisplayName', sprintf('KNSB | It=%.0fNs | Class %s', res_knsb.It, res_knsb.motor_class));
xlabel('Time [s]'); ylabel('Thrust [N]');
title('BATES grain: APCP vs KNSB', 'FontWeight','normal');
legend('Location','northeast');

subplot(1,2,2);
hold on; grid on; box on;
plot(res_apcp.t, res_apcp.Pc/1e6, 'b-', 'LineWidth', 2, 'DisplayName','APCP');
plot(res_knsb.t, res_knsb.Pc/1e6, 'r-', 'LineWidth', 2, 'DisplayName','KNSB');
xlabel('Time [s]'); ylabel('Chamber pressure [MPa]');
title('Chamber pressure comparison', 'FontWeight','normal');
legend('Location','northeast');
sgtitle('Phase 3 — Propellant Comparison on Fixed BATES Grain', ...
        'FontSize',12,'FontWeight','normal');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 4 — Trade studies + validation
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n--- PHASE 4: Trade Studies & Validation ---\n\n');

% ── Trade Study 1: Grain geometry comparison ──────────────────────────
fprintf('[Trade 1] Grain geometry comparison...\n');
T1 = trade_study('mode', 'grain_compare', ...
                 'propellant', 'APCP_Aerotech', ...
                 'Do', 0.075, 'Di', 0.025, 'L', 0.15, 'N', 4, ...
                 'At', At, 'epsilon', epsilon);

% ── Trade Study 2: Propellant comparison ─────────────────────────────
fprintf('\n[Trade 2] Propellant comparison on BATES grain...\n');
T2 = trade_study('mode', 'propellant_compare', ...
                 'grain_type', 'BATES', ...
                 'Do', 0.075, 'Di', 0.025, 'L', 0.15, 'N', 4, ...
                 'At', At, 'epsilon', epsilon);

% ── Trade Study 3: Parametric sweep — core diameter ──────────────────
fprintf('\n[Trade 3] BATES core diameter sweep...\n');
T3 = trade_study('mode', 'sweep', ...
                 'propellant', 'APCP_Aerotech', ...
                 'sweep_param', 'Di', ...
                 'sweep_range', [0.015, 0.060], ...
                 'n_sweep', 8, ...
                 'Do', 0.075, 'L', 0.15, 'N', 4, ...
                 'At', At, 'epsilon', epsilon);

% ── Validation: compare against built-in reference ───────────────────
fprintf('\n[Validation] Comparing APCP simulation vs. J350W reference...\n');
val = validation(res_apcp, 'builtin', 'plot', true);

%% ════════════════════════════════════════════════════════════════════════
%  FINAL SUMMARY TABLE
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n\n========== Final Simulation Summary ==========\n');
fprintf('Grain: BATES 75mm OD / 25mm core / 4 segments\n');
fprintf('Nozzle: Dt=%.1fmm  epsilon=%.0f  Sea level\n\n', nozzle_apcp.Dt*1e3, epsilon);

configs = {
    'APCP_Aerotech', res_apcp;
    'KNSB',          res_knsb;
};

fprintf('%-20s  %-6s  %-8s  %-8s  %-8s  %-6s  %-6s\n', ...
    'Propellant','Class','It[Ns]','Fpk[N]','Favg[N]','Pc[MPa]','Isp[s]');
fprintf('%s\n', repmat('-', 1, 72));
for i = 1:size(configs,1)
    name = configs{i,1};
    r    = configs{i,2};
    fprintf('%-20s  %-6s  %-8.1f  %-8.1f  %-8.1f  %-6.2f  %-6.1f\n', ...
        name, r.motor_class, r.It, r.F_max, r.F_avg, r.Pc_max/1e6, r.Isp);
end
fprintf('\nValidation (APCP vs J350W builtin): R²=%.4f  It_err=%.1f%%\n', ...
    val.r2, val.err_It);
fprintf('================================================\n\n');

fprintf('Phase 4 complete. All trade study figures generated.\n');
fprintf('Tip: Save figures with  saveas(gcf, ''filename.png'')\n');
fprintf('     or use report_gen.m for a full automated PDF export.\n');
