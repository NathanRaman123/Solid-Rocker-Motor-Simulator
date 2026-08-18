%% phase1_run.m — Phase 1 Demo: Grain Geometry Engine
%
%  Run this script to exercise all three grain types and compare their
%  burn profiles. No toolboxes required except Optimization Toolbox
%  (for fsolve in wagon_wheel.m — can substitute a manual bisection if needed).
%
%  Output: Figure with 2x2 panel comparing all three geometries.
%
%  Typical run time: < 1 second
%
%  After running, try modifying the parameters below and re-running
%  to build intuition for how geometry drives burn profile.
% ─────────────────────────────────────────────────────────────────────────

clear; clc; close all;

fprintf('=== Phase 1: Grain Geometry Engine ===\n\n');

%% ── 1. BATES grain ───────────────────────────────────────────────────
% Cylindrical core — simple, progressive, well-validated
% Think: AeroTech RMS reload hardware, most amateur HPR motors
Do_bates = 0.075;   % Outer diameter [m]  (75 mm — 3" motor)
Di_bates = 0.025;   % Core diameter  [m]  (25 mm)
L_bates  = 0.150;   % Segment length [m]  (150 mm)
N_bates  = 4;       % Number of segments

g_bates = BATES(Do_bates, Di_bates, L_bates, N_bates, 500);

%% ── 2. Star grain ────────────────────────────────────────────────────
% 6-point star — near-neutral burn, common in upper-stage motors
% epsilon controls how "filled in" the star is (higher = more circular)
Do_star   = 0.075;  % Match BATES OD for fair comparison
Ri_star   = 0.013;  % Tip-to-center radius [m]
epsilon   = 0.45;   % Star fraction (0.3=sharp star, 0.7=almost circular)
N_pts     = 6;      % Number of star points
L_star    = 0.20;   % Grain length [m]

g_star = star_grain(Do_star, Ri_star, epsilon, N_pts, L_star, 500);

%% ── 3. Wagon-wheel grain ─────────────────────────────────────────────
% Hub + 6 radial spokes — high initial Ab, fast pressurization
% Used when fast thrust rise (igniter charge, boost phase) is needed
Do_ww    = 0.080;   % Slightly larger motor [m]
Rc_ww    = 0.012;   % Hub radius [m]
N_spokes = 6;
spoke_w  = 0.006;   % Spoke width [m]
spoke_L  = 0.018;   % Spoke length from hub [m]
L_ww     = 0.20;    % Grain length [m]

g_ww = wagon_wheel(Do_ww, Rc_ww, N_spokes, spoke_w, spoke_L, L_ww, 500);

%% ── 4. Visualize all three ───────────────────────────────────────────
grain_viz(g_bates, g_star, g_ww);

%% ── 5. Print comparison table ────────────────────────────────────────
fprintf('\n%-20s %-15s %-12s %-12s %-10s\n', ...
    'Grain', 'Profile', 'Ab_init[cm2]', 'Ab_max[cm2]', 'Web[mm]');
fprintf('%s\n', repmat('-', 1, 72));

for g = {g_bates, g_star, g_ww}
    grain = g{1};
    fprintf('%-20s %-15s %-12.1f %-12.1f %-10.1f\n', ...
        grain.type, ...
        grain.profile_type, ...
        grain.Ab(1)   * 1e4, ...
        max(grain.Ab) * 1e4, ...
        grain.w_max   * 1e3);
end

fprintf('\n');
fprintf('Next step: feed Ab(t) and Vp(t) into Phase 2 propellant model\n');
fprintf('  >> propellant = propellant_db(''APCP_Aerotech'');\n');
fprintf('  >> [Pc, F, It] = ballistics_ode(g_bates, propellant, nozzle);\n');

%% ── 6. Sensitivity study: effect of BATES core diameter ──────────────
fprintf('\n[Sensitivity] BATES Ab_init vs. Di/Do ratio...\n');
Do_s   = 0.075;
Di_vec = linspace(0.015, 0.060, 10);
figure('Name','BATES Sensitivity — Di/Do vs. Burn Area', ...
       'Color',[1 1 1], 'Position',[200 200 700 400]);
hold on; grid on; box on;
for Di_s = Di_vec
    gs = BATES(Do_s, Di_s, 0.15, 4, 300);
    plot(gs.w*1e3, gs.Ab*1e4, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Di=%.0fmm (%.0f%%)', Di_s*1e3, Di_s/Do_s*100));
end
xlabel('Web regression [mm]');
ylabel('Burn area, A_b [cm²]');
title(sprintf('BATES sensitivity: core diameter Di (OD=%.0fmm fixed)', Do_s*1e3), ...
      'FontWeight','normal');
legend('Location','best','FontSize',8);
