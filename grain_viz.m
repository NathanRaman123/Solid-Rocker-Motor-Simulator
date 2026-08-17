function grain_viz(grain, varargin)
% GRAIN_VIZ  Visualization suite for solid rocket motor grain geometries
%
%   Generates a 2x2 figure panel:
%     (1) Cross-section at initial state
%     (2) Burn area Ab vs. web regression
%     (3) Propellant volume Vp vs. web regression
%     (4) Normalized burn area (Ab/Ab0) — useful for profile classification
%
%   Usage:
%     grain_viz(grain)                 % single grain
%     grain_viz(grain, grain2, grain3) % overlay multiple grains on plots 2-4
%
%   Compatible with outputs of BATES(), star_grain(), wagon_wheel()
%
%   Example (compare all three):
%     g1 = BATES(0.075, 0.025, 0.150, 4, 500);
%     g2 = star_grain(0.075, 0.015, 0.45, 6, 0.20, 500);
%     g3 = wagon_wheel(0.08, 0.012, 6, 0.006, 0.018, 0.20, 500);
%     grain_viz(g1, g2, g3);

    all_grains = [{grain}, varargin];
    n_grains   = numel(all_grains);
    colors     = lines(n_grains);

    fig = figure('Name', 'Grain Geometry Engine — Phase 1', ...
                 'NumberTitle', 'off', ...
                 'Color', [1 1 1], ...
                 'Position', [100 100 1100 820]);

    % ── Panel 1: Cross-section of primary grain ────────────────────────
    ax1 = subplot(2, 2, 1);
    draw_cross_section(ax1, grain);
    title(ax1, sprintf('%s grain — initial cross-section', grain.type), ...
          'FontWeight', 'normal');

    % ── Panel 2: Burn area vs. web ─────────────────────────────────────
    ax2 = subplot(2, 2, 2);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    for i = 1:n_grains
        g = all_grains{i};
        plot(ax2, g.w * 1e3, g.Ab * 1e4, ...
             'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', g.type);
    end
    xlabel(ax2, 'Web regression, w [mm]');
    ylabel(ax2, 'Burn area, A_b [cm²]');
    title(ax2, 'Burn area vs. web regression', 'FontWeight', 'normal');
    legend(ax2, 'Location', 'best');
    xlim(ax2, [0, max(cellfun(@(g) g.w_max, all_grains)) * 1e3 * 1.05]);

    % ── Panel 3: Propellant volume vs. web ────────────────────────────
    ax3 = subplot(2, 2, 3);
    hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
    for i = 1:n_grains
        g = all_grains{i};
        plot(ax3, g.w * 1e3, g.Vp * 1e6, ...
             'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', g.type);
    end
    xlabel(ax3, 'Web regression, w [mm]');
    ylabel(ax3, 'Propellant volume, V_p [cm³]');
    title(ax3, 'Propellant volume consumed vs. web', 'FontWeight', 'normal');
    legend(ax3, 'Location', 'best');

    % ── Panel 4: Normalized burn area ─────────────────────────────────
    ax4 = subplot(2, 2, 4);
    hold(ax4, 'on'); grid(ax4, 'on'); box(ax4, 'on');
    yline(ax4, 1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'Neutral (A_b = const)');
    for i = 1:n_grains
        g = all_grains{i};
        Ab0 = g.Ab(1);
        if Ab0 > 0
            w_norm  = g.w / g.w_max;
            Ab_norm = g.Ab / Ab0;
            plot(ax4, w_norm, Ab_norm, ...
                 'Color', colors(i,:), 'LineWidth', 2, ...
                 'DisplayName', sprintf('%s (%s)', g.type, g.profile_type));
        end
    end
    xlabel(ax4, 'Normalized web fraction, w/w_{max}');
    ylabel(ax4, 'Normalized burn area, A_b / A_{b,0}');
    title(ax4, 'Burn profile classification', 'FontWeight', 'normal');
    legend(ax4, 'Location', 'best', 'FontSize', 8);
    ylim(ax4, [0, max(2, ax4.YLim(2))]);

    % ── Global title ──────────────────────────────────────────────────
    sgtitle('Solid Rocket Motor — Grain Geometry Engine (Phase 1)', ...
            'FontSize', 13, 'FontWeight', 'normal');
end


% ─────────────────────────────────────────────────────────────────────────
function draw_cross_section(ax, grain)
% Draws the grain cross-section on axes ax.
    hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');
    R_OD = grain.Do / 2;
    theta = linspace(0, 2*pi, 360);

    % Case outer diameter
    fill(ax, R_OD*cos(theta), R_OD*sin(theta), [0.75 0.75 0.75], ...
         'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
         'DisplayName', 'Grain (propellant)');

    switch grain.type

        case 'BATES'
            % Core hole (white = void)
            Ri = grain.Di / 2;
            fill(ax, Ri*cos(theta), Ri*sin(theta), [1 1 1], ...
                 'EdgeColor', [0.2 0.2 0.8], 'LineWidth', 1.5, ...
                 'DisplayName', 'Port (void)');
            % Annotate web
            annotation_arc(ax, 0, Ri, R_OD, 'w_{radial}');

        case 'Star'
            cs = grain.cross_section;
            px = [cs.x, cs.x(1)];
            py = [cs.y, cs.y(1)];
            fill(ax, px, py, [1 1 1], ...
                 'EdgeColor', [0.2 0.2 0.8], 'LineWidth', 1.5, ...
                 'DisplayName', 'Port (void)');

        case 'WagonWheel'
            cs = grain.cross_section;
            % Hub
            fill(ax, cs.hub_x, cs.hub_y, [1 1 1], ...
                 'EdgeColor', [0.2 0.2 0.8], 'LineWidth', 1.5);
            % Spokes
            for j = 1:numel(cs.spk_x)
                fill(ax, cs.spk_x{j}, cs.spk_y{j}, [1 1 1], ...
                     'EdgeColor', [0.2 0.2 0.8], 'LineWidth', 1.2);
            end
    end

    % OD circle outline
    plot(ax, R_OD*cos(theta), R_OD*sin(theta), 'k-', 'LineWidth', 1.5);

    % Scale bar
    r_tick = R_OD * 1.15;
    text(ax, r_tick, 0, sprintf('OD = %.0f mm', grain.Do*1e3), ...
         'FontSize', 8, 'HorizontalAlignment', 'left', ...
         'Color', [0.3 0.3 0.3]);

    axis(ax, [-R_OD*1.3, R_OD*1.3, -R_OD*1.3, R_OD*1.3]);
    xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]');
end


function annotation_arc(ax, ~, r_inner, r_outer, label)
% Small annotation showing web thickness.
    mid_r = (r_inner + r_outer) / 2;
    text(ax, mid_r * 0.7, mid_r * 0.7, label, ...
         'FontSize', 8, 'Color', [0.2 0.5 0.2], ...
         'HorizontalAlignment', 'center');
end
