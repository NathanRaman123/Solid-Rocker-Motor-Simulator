%% SRM_SIMULATOR.m — Solid Rocket Motor Simulator (All Phases)
%
%  Single-file master runner. Executes the full simulation pipeline:
%
%    Phase 1 — Grain geometry engine (BATES, Star, Wagon-Wheel)
%    Phase 2 — Propellant thermochemistry + nozzle model
%    Phase 3 — Internal ballistics ODE (ode45)
%    Phase 4 — Trade studies, parametric sweeps, validation
%
%  Prerequisites: All supporting .m files in the same folder as this script.
%  Toolboxes:     Optimization Toolbox (fsolve). See bottom for no-toolbox fix.
%
%  Run time: ~15–45 seconds depending on machine.
%  Figures:  8 figures produced, each in its own named window.
%
%  Usage:
%    >> SRM_SIMULATOR
% ─────────────────────────────────────────────────────────────────────────

clear; clc; close all;

fprintf('╔══════════════════════════════════════════════════════╗\n');
fprintf('║   Solid Rocket Motor Simulator — Full Pipeline       ║\n');
fprintf('║                    Purdue Aerospace                  ║\n');
fprintf('╚══════════════════════════════════════════════════════╝\n\n');

%% ════════════════════════════════════════════════════════════════════════
%  SHARED DESIGN PARAMETERS
%  (Change these to explore different motor configurations)
%% ════════════════════════════════════════════════════════════════════════

% --- Grain (BATES reference design) ---
Do      = 0.075;   % Outer diameter [m]          — 75 mm (3" motor)
Di      = 0.025;   % Core diameter  [m]          — 25 mm
L_seg   = 0.150;   % Segment length [m]
N_seg   = 4;       % Number of segments

% --- Nozzle ---
Dt      = 0.015;   % Throat diameter [m]         — 15 mm
At      = pi*(Dt/2)^2;
epsilon = 8;       % Expansion ratio Ae/At

% --- Environment ---
P_amb   = 101325;  % Ambient pressure [Pa]       — sea level

n_pts   = 500;     % Discretization steps for grain geometry

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 1 — GRAIN GEOMETRY ENGINE
%% ════════════════════════════════════════════════════════════════════════

fprintf('▶  PHASE 1: Grain Geometry Engine\n');
fprintf('   ─────────────────────────────────────────────────────\n');

% ── Build all three grain types ───────────────────────────────────────

g_bates = BATES(Do, Di, L_seg, N_seg, n_pts);

Ri_star   = Di/2 * 0.85;          % Star tip radius — slightly inside core
g_star    = star_grain(Do, Ri_star, 0.45, 6, L_seg*N_seg, n_pts);

Rc_ww     = Di/2 * 0.70;
spk_w     = Do * 0.06;
spk_L     = Do * 0.16;
g_ww      = wagon_wheel(Do, Rc_ww, 6, spk_w, spk_L, L_seg*N_seg, n_pts);

% ── FIGURE 1: Grain cross-sections + burn profiles ────────────────────

fprintf('\n   Generating Figure 1: Grain geometry comparison...\n');

fig1 = figure('Name','[P1] Grain Geometry Engine','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[30 430 1300 620]);

all_grains = {g_bates, g_star, g_ww};
grain_colors = {[0.18 0.44 0.71],[0.88 0.30 0.20],[0.20 0.63 0.17]};

% Row 1: cross-sections
for k = 1:3
    subplot(2,3,k);
    draw_cross_section_inline(all_grains{k}, grain_colors{k});
end

% Row 2: burn profiles overlay on shared axes
ax_ab = subplot(2,3,4);
hold(ax_ab,'on'); grid(ax_ab,'on'); box(ax_ab,'on');
for k = 1:3
    g = all_grains{k};
    plot(ax_ab, g.w*1e3, g.Ab*1e4, 'Color', grain_colors{k}, ...
         'LineWidth', 2.2, 'DisplayName', g.type);
end
xlabel(ax_ab,'Web regression [mm]'); ylabel(ax_ab,'Burn area A_b [cm²]');
title(ax_ab,'Burn area vs. web','FontWeight','normal');
legend(ax_ab,'Location','best');

ax_vp = subplot(2,3,5);
hold(ax_vp,'on'); grid(ax_vp,'on'); box(ax_vp,'on');
for k = 1:3
    g = all_grains{k};
    plot(ax_vp, g.w*1e3, g.Vp*1e6, 'Color', grain_colors{k}, ...
         'LineWidth', 2.2, 'DisplayName', g.type);
end
xlabel(ax_vp,'Web regression [mm]'); ylabel(ax_vp,'Propellant volume [cm³]');
title(ax_vp,'Propellant volume vs. web','FontWeight','normal');
legend(ax_vp,'Location','best');

ax_norm = subplot(2,3,6);
hold(ax_norm,'on'); grid(ax_norm,'on'); box(ax_norm,'on');
yline(ax_norm, 1.0, 'k--', 'LineWidth', 1.2, 'DisplayName','Neutral');
for k = 1:3
    g = all_grains{k};
    Ab0 = g.Ab(1);
    if Ab0 > 0
        plot(ax_norm, g.w/g.w_max, g.Ab/Ab0, ...
             'Color', grain_colors{k}, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('%s (%s)', g.type, g.profile_type));
    end
end
xlabel(ax_norm,'Normalized web w/w_{max}');
ylabel(ax_norm,'Normalized A_b / A_{b,0}');
title(ax_norm,'Burn profile classification','FontWeight','normal');
legend(ax_norm,'Location','best','FontSize',8);

sgtitle(fig1,'Phase 1 — Grain Geometry Engine','FontSize',14,'FontWeight','normal');

% ── FIGURE 2: BATES sensitivity — core diameter sweep ─────────────────

fprintf('   Generating Figure 2: BATES core diameter sensitivity...\n');

fig2 = figure('Name','[P1] BATES Sensitivity','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[30 30 750 400]);
hold on; grid on; box on;
Di_vec = linspace(0.015, 0.060, 10);
cmap   = cool(numel(Di_vec));
for idx = 1:numel(Di_vec)
    gs = BATES(Do, Di_vec(idx), L_seg, N_seg, 300);
    plot(gs.w*1e3, gs.Ab*1e4, 'Color', cmap(idx,:), 'LineWidth', 1.6, ...
         'DisplayName', sprintf('Di=%.0fmm (%.0f%%)', Di_vec(idx)*1e3, Di_vec(idx)/Do*100));
end
xlabel('Web regression [mm]'); ylabel('Burn area A_b [cm²]');
title(sprintf('BATES sensitivity: core diameter  (OD = %.0f mm fixed)', Do*1e3), ...
      'FontWeight','normal');
legend('Location','best','FontSize',7.5);
colorbar('Ticks',[0 1],'TickLabels',{'Di small','Di large'});

% ── Print Phase 1 table ───────────────────────────────────────────────

fprintf('\n   %-20s %-18s %-13s %-13s %-10s\n', ...
    'Grain','Profile','Ab_init[cm²]','Ab_max[cm²]','Web[mm]');
fprintf('   %s\n', repmat('-',1,76));
for k = 1:3
    g = all_grains{k};
    fprintf('   %-20s %-18s %-13.1f %-13.1f %-10.1f\n', ...
        g.type, g.profile_type, g.Ab(1)*1e4, max(g.Ab)*1e4, g.w_max*1e3);
end
fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 2 — PROPELLANT THERMOCHEMISTRY + NOZZLE MODEL
%% ════════════════════════════════════════════════════════════════════════

fprintf('▶  PHASE 2: Propellant & Nozzle Model\n');
fprintf('   ─────────────────────────────────────────────────────\n');

prop_apcp = propellant_db('APCP_Aerotech');
prop_knsb = propellant_db('KNSB');

nozzle_apcp = nozzle_model(prop_apcp, At, epsilon, P_amb);
nozzle_knsb = nozzle_model(prop_knsb, At, epsilon, P_amb);

[~, rc_apcp] = burn_rate(prop_apcp, []);
[~, rc_knsb] = burn_rate(prop_knsb, []);

% ── FIGURE 3: Propellant characterization ─────────────────────────────

fprintf('\n   Generating Figure 3: Propellant characterization...\n');

fig3 = figure('Name','[P2] Propellant & Nozzle','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[380 430 1200 500]);

% Panel 1: Burn rate vs pressure with temperature bands
ax31 = subplot(1,3,1);
hold(ax31,'on'); grid(ax31,'on'); box(ax31,'on');
fill(ax31,[rc_apcp.P, fliplr(rc_apcp.P)]/1e6, ...
     [rc_apcp.r_cold, fliplr(rc_apcp.r_hot)]*1e3, ...
     [0.78 0.86 1.0],'EdgeColor','none','FaceAlpha',0.5,'DisplayName','APCP ±30K');
fill(ax31,[rc_knsb.P, fliplr(rc_knsb.P)]/1e6, ...
     [rc_knsb.r_cold, fliplr(rc_knsb.r_hot)]*1e3, ...
     [1.0 0.82 0.78],'EdgeColor','none','FaceAlpha',0.5,'DisplayName','KNSB ±30K');
plot(ax31, rc_apcp.P/1e6, rc_apcp.r*1e3, 'b-', 'LineWidth',2.5, 'DisplayName','APCP Aerotech');
plot(ax31, rc_knsb.P/1e6, rc_knsb.r*1e3, 'r-', 'LineWidth',2.5, 'DisplayName','KNSB');
xlabel(ax31,'Chamber pressure [MPa]'); ylabel(ax31,'Burn rate [mm/s]');
title(ax31,"Saint-Robert's law  r = a·P^n",'FontWeight','normal');
legend(ax31,'Location','northwest','FontSize',8); xlim(ax31,[0 12]);

% Panel 2: Log-log (slope = n)
ax32 = subplot(1,3,2);
hold(ax32,'on'); grid(ax32,'on'); box(ax32,'on');
loglog(ax32, rc_apcp.P/1e6, rc_apcp.r*1e3, 'b-', 'LineWidth',2.5, ...
       'DisplayName',sprintf('APCP (n=%.3f)',prop_apcp.n));
loglog(ax32, rc_knsb.P/1e6, rc_knsb.r*1e3, 'r-', 'LineWidth',2.5, ...
       'DisplayName',sprintf('KNSB (n=%.3f)',prop_knsb.n));
xlabel(ax32,'Pressure [MPa]'); ylabel(ax32,'Burn rate [mm/s]');
title(ax32,'Log-log  (slope = pressure exponent n)','FontWeight','normal');
legend(ax32,'Location','northwest','FontSize',8);

% Panel 3: Thrust coefficient vs Pc
ax33 = subplot(1,3,3);
hold(ax33,'on'); grid(ax33,'on'); box(ax33,'on');
Pc_range = linspace(0.5e6, 8e6, 200);
Cf_apcp  = arrayfun(@(P) nozzle_apcp.Cf(P), Pc_range);
Cf_knsb  = arrayfun(@(P) nozzle_knsb.Cf(P), Pc_range);
plot(ax33, Pc_range/1e6, Cf_apcp, 'b-', 'LineWidth',2.5, ...
     'DisplayName',sprintf('APCP  ε=%.0f',epsilon));
plot(ax33, Pc_range/1e6, Cf_knsb, 'r-', 'LineWidth',2.5, ...
     'DisplayName',sprintf('KNSB  ε=%.0f',epsilon));
xlabel(ax33,'Chamber pressure [MPa]'); ylabel(ax33,'Thrust coefficient C_f');
title(ax33,sprintf('Nozzle C_f  (D_t = %.0f mm)', Dt*1e3),'FontWeight','normal');
legend(ax33,'Location','southeast','FontSize',8);

sgtitle(fig3,'Phase 2 — Propellant Thermochemistry & Nozzle','FontSize',14,'FontWeight','normal');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 3 — INTERNAL BALLISTICS ODE
%% ════════════════════════════════════════════════════════════════════════

fprintf('\n▶  PHASE 3: Internal Ballistics ODE Solver\n');
fprintf('   ─────────────────────────────────────────────────────\n');

res_apcp = ballistics_ode(g_bates, prop_apcp, nozzle_apcp);
res_knsb = ballistics_ode(g_bates, prop_knsb, nozzle_knsb);

% Also run star and wagon-wheel with APCP for grain comparison in Phase 4
res_star = ballistics_ode(g_star, prop_apcp, nozzle_apcp);
res_ww   = ballistics_ode(g_ww,   prop_apcp, nozzle_apcp);

% ── FIGURE 4: Full ballistics dashboard (BATES/APCP reference design) ─

fprintf('\n   Generating Figure 4: Full ballistics dashboard...\n');

fig4 = figure('Name','[P3] Ballistics Dashboard — BATES/APCP','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[30 30 1250 720]);

t  = res_apcp.t;
tb = res_apcp.tb;

% Thrust
ax41 = subplot(2,3,1);
plot(ax41, t, res_apcp.F, 'b-', 'LineWidth',2.5);
hold(ax41,'on'); grid(ax41,'on'); box(ax41,'on');
xline(ax41, tb,'r--','Burnout','LineWidth',1.2,'LabelOrientation','horizontal');
area(ax41, t, res_apcp.F, 'FaceColor',[0.75 0.87 1.0],'EdgeColor','none','FaceAlpha',0.4);
xlabel(ax41,'Time [s]'); ylabel(ax41,'Thrust [N]');
title(ax41,sprintf('Thrust  (Class %s | I_t = %.0f N·s)', ...
      res_apcp.motor_class, res_apcp.It),'FontWeight','normal');

% Chamber pressure
ax42 = subplot(2,3,2);
plot(ax42, t, res_apcp.Pc/1e6, 'r-', 'LineWidth',2.5);
hold(ax42,'on'); grid(ax42,'on'); box(ax42,'on');
xline(ax42, tb,'r--','LineWidth',1.2);
yline(ax42, res_apcp.Pc_max/1e6,'k:','LineWidth',1.2);
xlabel(ax42,'Time [s]'); ylabel(ax42,'P_c [MPa]');
title(ax42,sprintf('Chamber pressure  (peak %.2f MPa)', ...
      res_apcp.Pc_max/1e6),'FontWeight','normal');
text(ax42, 0.97, 0.92, sprintf('P_{max} = %.2f MPa', res_apcp.Pc_max/1e6), ...
     'Units','normalized','HorizontalAlignment','right','FontSize',8);

% Web regression over time
ax43 = subplot(2,3,3);
plot(ax43, t, res_apcp.w*1e3,'g-','LineWidth',2.5);
hold(ax43,'on'); grid(ax43,'on'); box(ax43,'on');
yline(ax43, g_bates.w_max*1e3,'k--','w_{max}','LineWidth',1.2);
xlabel(ax43,'Time [s]'); ylabel(ax43,'Web regression [mm]');
title(ax43,'Web regression w(t)','FontWeight','normal');

% Burn rate vs pressure (phase portrait colored by time)
ax44 = subplot(2,3,4);
scatter(ax44, res_apcp.Pc/1e6, res_apcp.r*1e3, 6, t, 'filled');
hold(ax44,'on'); grid(ax44,'on'); box(ax44,'on');
P_line = linspace(min(res_apcp.Pc), max(res_apcp.Pc), 200);
plot(ax44, P_line/1e6, prop_apcp.a.*P_line.^prop_apcp.n*1e3, ...
     'k-','LineWidth',1.5,'DisplayName','r = aP^n');
cb = colorbar(ax44); cb.Label.String = 'Time [s]';
xlabel(ax44,'Pressure [MPa]'); ylabel(ax44,'Burn rate [mm/s]');
title(ax44,sprintf("Saint-Robert's  n = %.3f", prop_apcp.n),'FontWeight','normal');

% Burn area vs web
ax45 = subplot(2,3,5);
plot(ax45, res_apcp.w*1e3, res_apcp.Ab*1e4, 'm-', 'LineWidth',2.5);
grid(ax45,'on'); box(ax45,'on');
xlabel(ax45,'Web [mm]'); ylabel(ax45,'A_b [cm²]');
title(ax45,'Burn area A_b(w)','FontWeight','normal');

% Performance summary text panel
ax46 = subplot(2,3,6);
axis(ax46,'off');
lines_txt = {
    sprintf('┌─ Motor Performance ──────────┐');
    sprintf('  Class         :  %s',          res_apcp.motor_class);
    sprintf('  Total impulse :  %.1f N·s',    res_apcp.It);
    sprintf('  Burn duration :  %.3f s',      res_apcp.tb);
    sprintf('  Peak thrust   :  %.1f N',      res_apcp.F_max);
    sprintf('  Avg thrust    :  %.1f N',      res_apcp.F_avg);
    sprintf('  Peak Pc       :  %.2f MPa',    res_apcp.Pc_max/1e6);
    sprintf('  Delivered Isp :  %.1f s',      res_apcp.Isp);
    sprintf('  Prop mass     :  %.3f kg',     res_apcp.mp);
    sprintf('├─ Design Inputs ──────────────┤');
    sprintf('  Grain         :  %s', g_bates.type);
    sprintf('  Propellant    :  %s', prop_apcp.name);
    sprintf('  OD / Core     :  %.0f / %.0f mm', Do*1e3, Di*1e3);
    sprintf('  Throat Dt     :  %.1f mm', Dt*1e3);
    sprintf('  Expansion ε   :  %.0f', epsilon);
    sprintf('└──────────────────────────────┘');
};
for i = 1:numel(lines_txt)
    text(ax46, 0.03, 1.0-(i-1)*0.063, lines_txt{i}, ...
         'Units','normalized','FontSize',8.5,'FontName','Monospaced', ...
         'VerticalAlignment','top');
end

sgtitle(fig4,'Phase 3 — Internal Ballistics  (BATES / APCP Aerotech)', ...
        'FontSize',14,'FontWeight','normal');

% ── FIGURE 5: APCP vs KNSB head-to-head ──────────────────────────────

fprintf('   Generating Figure 5: APCP vs KNSB comparison...\n');

fig5 = figure('Name','[P3] Propellant Head-to-Head','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[380 30 1000 430]);

subplot(1,2,1);
hold on; grid on; box on;
plot(res_apcp.t, res_apcp.F, 'b-', 'LineWidth',2.5, ...
     'DisplayName', sprintf('APCP  |  It=%.0f Ns  Class %s', res_apcp.It, res_apcp.motor_class));
plot(res_knsb.t, res_knsb.F, 'r-', 'LineWidth',2.5, ...
     'DisplayName', sprintf('KNSB  |  It=%.0f Ns  Class %s', res_knsb.It, res_knsb.motor_class));
xlabel('Time [s]'); ylabel('Thrust [N]');
title('Thrust — same BATES grain, two propellants','FontWeight','normal');
legend('Location','northeast','FontSize',9);

subplot(1,2,2);
hold on; grid on; box on;
plot(res_apcp.t, res_apcp.Pc/1e6, 'b-', 'LineWidth',2.5, 'DisplayName','APCP');
plot(res_knsb.t, res_knsb.Pc/1e6, 'r-', 'LineWidth',2.5, 'DisplayName','KNSB');
xlabel('Time [s]'); ylabel('Chamber pressure [MPa]');
title('Chamber pressure comparison','FontWeight','normal');
legend('Location','northeast');

sgtitle(fig5,'Phase 3 — Propellant Comparison on Fixed BATES Grain', ...
        'FontSize',13,'FontWeight','normal');

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 4 — TRADE STUDIES + VALIDATION
%% ════════════════════════════════════════════════════════════════════════

fprintf('\n▶  PHASE 4: Trade Studies & Validation\n');
fprintf('   ─────────────────────────────────────────────────────\n');

% ── FIGURE 6: Grain geometry trade study ─────────────────────────────

fprintf('\n   Generating Figure 6: Grain geometry trade study...\n');

all_res    = {res_apcp, res_star, res_ww};
all_labels = {'BATES','Star (6-pt)','WagonWheel'};
trade_colors = {[0.18 0.44 0.71],[0.88 0.30 0.20],[0.20 0.63 0.17]};

fig6 = figure('Name','[P4] Grain Trade Study','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[30 30 1250 700]);

% Thrust overlay
ax61 = subplot(2,3,1);
hold(ax61,'on'); grid(ax61,'on'); box(ax61,'on');
for k = 1:3
    r = all_res{k};
    plot(ax61, r.t, r.F, 'Color',trade_colors{k},'LineWidth',2.2, ...
         'DisplayName', sprintf('%s (Class %s)', all_labels{k}, r.motor_class));
end
xlabel(ax61,'Time [s]'); ylabel(ax61,'Thrust [N]');
title(ax61,'Thrust curves — grain comparison','FontWeight','normal');
legend(ax61,'Location','northeast','FontSize',8);

% Pressure overlay
ax62 = subplot(2,3,2);
hold(ax62,'on'); grid(ax62,'on'); box(ax62,'on');
for k = 1:3
    r = all_res{k};
    plot(ax62, r.t, r.Pc/1e6, 'Color',trade_colors{k},'LineWidth',2.2, ...
         'DisplayName', all_labels{k});
end
xlabel(ax62,'Time [s]'); ylabel(ax62,'P_c [MPa]');
title(ax62,'Chamber pressure — grain comparison','FontWeight','normal');
legend(ax62,'Location','northeast','FontSize',8);

% Total impulse bar
ax63 = subplot(2,3,3);
It_vals = cellfun(@(r) r.It, all_res);
bh = bar(ax63, It_vals, 0.55);
bh.FaceColor = 'flat';
for k=1:3, bh.CData(k,:) = trade_colors{k}; end
set(ax63,'XTick',1:3,'XTickLabel',all_labels,'XTickLabelRotation',12);
ylabel(ax63,'Total impulse [N·s]');
title(ax63,'Total impulse','FontWeight','normal');
grid(ax63,'on'); box(ax63,'on');
for k=1:3
    text(ax63, k, It_vals(k)+5, sprintf('%.0f',It_vals(k)), ...
         'HorizontalAlignment','center','FontSize',9);
end

% Peak pressure bar
ax64 = subplot(2,3,4);
Pc_vals = cellfun(@(r) r.Pc_max/1e6, all_res);
bh2 = bar(ax64, Pc_vals, 0.55);
bh2.FaceColor = 'flat';
for k=1:3, bh2.CData(k,:) = trade_colors{k}; end
set(ax64,'XTick',1:3,'XTickLabel',all_labels,'XTickLabelRotation',12);
ylabel(ax64,'Peak P_c [MPa]');
title(ax64,'Peak chamber pressure','FontWeight','normal');
grid(ax64,'on'); box(ax64,'on');

% Delivered Isp bar
ax65 = subplot(2,3,5);
Isp_vals = cellfun(@(r) r.Isp, all_res);
bh3 = bar(ax65, Isp_vals, 0.55);
bh3.FaceColor = 'flat';
for k=1:3, bh3.CData(k,:) = trade_colors{k}; end
set(ax65,'XTick',1:3,'XTickLabel',all_labels,'XTickLabelRotation',12);
ylabel(ax65,'Delivered I_{sp} [s]');
title(ax65,'Specific impulse','FontWeight','normal');
grid(ax65,'on'); box(ax65,'on');

% Burn duration bar
ax66 = subplot(2,3,6);
tb_vals = cellfun(@(r) r.tb, all_res);
bh4 = bar(ax66, tb_vals, 0.55);
bh4.FaceColor = 'flat';
for k=1:3, bh4.CData(k,:) = trade_colors{k}; end
set(ax66,'XTick',1:3,'XTickLabel',all_labels,'XTickLabelRotation',12);
ylabel(ax66,'Burn duration [s]');
title(ax66,'Burn duration','FontWeight','normal');
grid(ax66,'on'); box(ax66,'on');

sgtitle(fig6,'Phase 4 — Grain Geometry Trade Study  (APCP Aerotech, fixed nozzle)', ...
        'FontSize',14,'FontWeight','normal');

% ── FIGURE 7: Parametric sweep — BATES core diameter ─────────────────

fprintf('   Generating Figure 7: Parametric sweep (core diameter)...\n');

Di_sweep  = linspace(0.015, 0.058, 9);
cmap_sw   = parula(numel(Di_sweep));
sweep_res = cell(1,numel(Di_sweep));

for idx = 1:numel(Di_sweep)
    gs = BATES(Do, Di_sweep(idx), L_seg, N_seg, 400);
    sweep_res{idx} = ballistics_ode(gs, prop_apcp, nozzle_apcp);
    sweep_res{idx}.Di_val = Di_sweep(idx);
end

fig7 = figure('Name','[P4] Parametric Sweep — Core Diameter','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[30 30 1200 500]);

ax71 = subplot(1,3,1);
hold(ax71,'on'); grid(ax71,'on'); box(ax71,'on');
for idx = 1:numel(Di_sweep)
    r = sweep_res{idx};
    plot(ax71, r.t, r.F, 'Color',cmap_sw(idx,:),'LineWidth',1.8, ...
         'DisplayName', sprintf('Di=%.0fmm', Di_sweep(idx)*1e3));
end
xlabel(ax71,'Time [s]'); ylabel(ax71,'Thrust [N]');
title(ax71,'Thrust vs. core diameter','FontWeight','normal');
legend(ax71,'Location','best','FontSize',7);
cb71 = colorbar(ax71,'Location','eastoutside');
cb71.Label.String = 'Di [mm]';
colormap(ax71, parula);
clim(ax71, [Di_sweep(1) Di_sweep(end)]*1e3);

ax72 = subplot(1,3,2);
hold(ax72,'on'); grid(ax72,'on'); box(ax72,'on');
It_sw  = cellfun(@(r) r.It,  sweep_res);
Pc_sw  = cellfun(@(r) r.Pc_max/1e6, sweep_res);
yyaxis(ax72,'left');
plot(ax72, Di_sweep*1e3, It_sw, 'b-o', 'LineWidth',2, 'MarkerSize',6);
ylabel(ax72,'Total impulse [N·s]');
yyaxis(ax72,'right');
plot(ax72, Di_sweep*1e3, Pc_sw, 'r-s', 'LineWidth',2, 'MarkerSize',6);
ylabel(ax72,'Peak P_c [MPa]');
xlabel(ax72,'Core diameter Di [mm]');
title(ax72,'Impulse & peak pressure vs. Di','FontWeight','normal');
legend(ax72,{'Total impulse','Peak Pc'},'Location','best');

ax73 = subplot(1,3,3);
hold(ax73,'on'); grid(ax73,'on'); box(ax73,'on');
F_avg_sw = cellfun(@(r) r.F_avg, sweep_res);
tb_sw    = cellfun(@(r) r.tb,   sweep_res);
yyaxis(ax73,'left');
plot(ax73, Di_sweep*1e3, F_avg_sw,'b-o','LineWidth',2,'MarkerSize',6);
ylabel(ax73,'Avg thrust [N]');
yyaxis(ax73,'right');
plot(ax73, Di_sweep*1e3, tb_sw,'r-s','LineWidth',2,'MarkerSize',6);
ylabel(ax73,'Burn duration [s]');
xlabel(ax73,'Core diameter Di [mm]');
title(ax73,'Avg thrust & burn duration vs. Di','FontWeight','normal');
legend(ax73,{'Avg thrust','Burn duration'},'Location','best');

sgtitle(fig7,'Phase 4 — Parametric Sweep: BATES Core Diameter  (OD fixed)', ...
        'FontSize',14,'FontWeight','normal');

% ── FIGURE 8: Validation ─────────────────────────────────────────────

fprintf('   Generating Figure 8: Validation against reference data...\n');

val = validation(res_apcp, 'builtin', 'plot', false);

fig8 = figure('Name','[P4] Validation','NumberTitle','off', ...
              'Color',[1 1 1],'Position',[380 30 1050 480]);

ax81 = subplot(1,2,1);
hold(ax81,'on'); grid(ax81,'on'); box(ax81,'on');
plot(ax81, val.t_ref, val.F_ref, 'k-','LineWidth',3.0, ...
     'DisplayName', sprintf('Reference: %s', val.motor_name));
plot(ax81, val.t_ref, val.F_sim, 'b--','LineWidth',2.2, ...
     'DisplayName', sprintf('Simulation  (R²=%.4f)', val.r2));
fill(ax81, [val.t_ref; flipud(val.t_ref)], ...
     [val.F_ref; flipud(val.F_sim)], ...
     [1 0.80 0.78],'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
xlabel(ax81,'Time [s]'); ylabel(ax81,'Thrust [N]');
title(ax81,'Thrust: simulation vs. reference','FontWeight','normal');
legend(ax81,'Location','northeast','FontSize',9);
text(ax81, 0.97, 0.45, ...
     sprintf('I_t err: %+.1f%%\nF_{pk} err: %+.1f%%\nRMSE: %.1f N\nR²: %.4f', ...
             val.err_It, val.err_Fpeak, val.rmse, val.r2), ...
     'Units','normalized','HorizontalAlignment','right','FontSize',9, ...
     'BackgroundColor',[0.96 0.96 0.96],'EdgeColor',[0.8 0.8 0.8]);

ax82 = subplot(1,2,2);
hold(ax82,'on'); grid(ax82,'on'); box(ax82,'on');
resid = val.F_sim - val.F_ref;
plot(ax82, val.t_ref, resid,'r-','LineWidth',1.8);
yline(ax82,  0,        'k-', 'LineWidth',1.2);
yline(ax82,  val.rmse, 'b--','±RMSE','LineWidth',1.2,'LabelHorizontalAlignment','left');
yline(ax82, -val.rmse, 'b--','LineWidth',1.2);
fill(ax82, [val.t_ref; flipud(val.t_ref)], ...
     [val.F_ref*0+val.rmse; flipud(val.F_ref*0-val.rmse)], ...
     [0.80 0.88 1.0],'EdgeColor','none','FaceAlpha',0.3);
xlabel(ax82,'Time [s]'); ylabel(ax82,'Residual F_{sim} - F_{ref} [N]');
title(ax82,'Thrust residuals','FontWeight','normal');

sgtitle(fig8,sprintf('Phase 4 — Validation  (%s)', val.motor_name), ...
        'FontSize',14,'FontWeight','normal');

%% ════════════════════════════════════════════════════════════════════════
%  FINAL SUMMARY
%% ════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                   FINAL SIMULATION SUMMARY                      ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Design: BATES %.0f/%.0f mm  |  Nozzle Dt=%.0fmm  ε=%.0f  |  Sea level   ║\n', ...
        Do*1e3, Di*1e3, Dt*1e3, epsilon);
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  %-14s  %-6s  %-9s  %-8s  %-8s  %-6s  %-5s  ║\n', ...
        'Propellant','Class','It [N·s]','Fpk [N]','Favg [N]','Pc[MPa]','Isp');
fprintf('║  %s  ║\n', repmat('─',1,62));

configs = {'APCP_Aerotech', res_apcp; 'KNSB', res_knsb};
for i = 1:size(configs,1)
    r = configs{i,2};
    fprintf('║  %-14s  %-6s  %-9.1f  %-8.1f  %-8.1f  %-6.2f  %-5.1f  ║\n', ...
        configs{i,1}, r.motor_class, r.It, r.F_max, r.F_avg, r.Pc_max/1e6, r.Isp);
end

fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Grain trade study (APCP):                                       ║\n');
grain_names_ts = {'BATES','Star','WagonWheel'};
for k = 1:3
    r = all_res{k};
    fprintf('║    %-12s  Class %-2s  It=%-7.0fNs  Fpk=%-6.0fN  tb=%.2fs        ║\n', ...
        grain_names_ts{k}, r.motor_class, r.It, r.F_max, r.tb);
end

fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Validation vs. J350W builtin reference:                         ║\n');
fprintf('║    R² = %.4f   |   It error = %+.1f%%   |   RMSE = %.1f N           ║\n', ...
        val.r2, val.err_It, val.rmse);
fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

fprintf('  8 figures generated:\n');
fprintf('    Fig 1  [P1] Grain geometry + burn profiles\n');
fprintf('    Fig 2  [P1] BATES core diameter sensitivity\n');
fprintf('    Fig 3  [P2] Propellant characterization + nozzle Cf\n');
fprintf('    Fig 4  [P3] Full ballistics dashboard (BATES/APCP)\n');
fprintf('    Fig 5  [P3] APCP vs KNSB head-to-head\n');
fprintf('    Fig 6  [P4] Grain geometry trade study\n');
fprintf('    Fig 7  [P4] Parametric sweep — core diameter\n');
fprintf('    Fig 8  [P4] Validation vs. reference thrust curve\n\n');

fprintf('  To save all figures:  savefigs_all\n');
fprintf('  To generate report:   report_gen\n\n');

%% ════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER — cross-section renderer (inline, no external dependency)
%% ════════════════════════════════════════════════════════════════════════

function draw_cross_section_inline(grain, clr)
    hold on; axis equal; grid on; box on;
    R_OD   = grain.Do / 2;
    theta  = linspace(0, 2*pi, 360);

    % Propellant body
    fill(R_OD*cos(theta), R_OD*sin(theta), clr, ...
         'FaceAlpha', 0.55, 'EdgeColor', clr*0.6, 'LineWidth', 1.5);

    switch grain.type

        case 'BATES'
            Ri = grain.Di / 2;
            fill(Ri*cos(theta), Ri*sin(theta), [1 1 1], ...
                 'EdgeColor',[0.25 0.25 0.75],'LineWidth',1.8);
            text(0, 0, 'core', 'HorizontalAlignment','center', ...
                 'FontSize',8,'Color',[0.3 0.3 0.3]);

        case 'Star'
            cs = grain.cross_section;
            px = [cs.x, cs.x(1)];
            py = [cs.y, cs.y(1)];
            fill(px, py, [1 1 1], ...
                 'EdgeColor',[0.25 0.25 0.75],'LineWidth',1.8);

        case 'WagonWheel'
            cs = grain.cross_section;
            fill(cs.hub_x, cs.hub_y, [1 1 1], ...
                 'EdgeColor',[0.25 0.25 0.75],'LineWidth',1.8);
            for j = 1:numel(cs.spk_x)
                fill(cs.spk_x{j}, cs.spk_y{j}, [1 1 1], ...
                     'EdgeColor',[0.25 0.25 0.75],'LineWidth',1.2);
            end
    end

    plot(R_OD*cos(theta), R_OD*sin(theta), 'k-', 'LineWidth', 2.0);
    axis([-R_OD*1.25, R_OD*1.25, -R_OD*1.25, R_OD*1.25]);
    xlabel('x [m]'); ylabel('y [m]');
    title(sprintf('%s — OD %.0fmm', grain.type, grain.Do*1e3), ...
          'FontWeight','normal');
end
