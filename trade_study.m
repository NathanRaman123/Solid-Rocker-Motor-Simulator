function results_table = trade_study(varargin)
% TRADE_STUDY  Parametric trade study across grain geometries and propellants
%
%   Sweeps user-defined parameter ranges and compares motor performance
%   across grain types and/or propellants. Generates a results table
%   and multi-panel comparison plots.
%
%   Usage (two modes):
%
%   MODE 1 — Grain geometry comparison (fixed propellant, vary grain):
%     results_table = trade_study('mode', 'grain_compare', ...
%                                 'propellant', 'APCP_Aerotech', ...
%                                 'At', pi*(0.015/2)^2, ...
%                                 'epsilon', 8, ...
%                                 'P_amb', 101325)
%
%   MODE 2 — Propellant comparison (fixed geometry, vary propellant):
%     results_table = trade_study('mode', 'propellant_compare', ...
%                                 'grain_type', 'BATES', ...
%                                 'Do', 0.075, 'Di', 0.025, ...
%                                 'L', 0.15, 'N', 4, ...
%                                 'At', pi*(0.015/2)^2, ...
%                                 'epsilon', 8)
%
%   MODE 3 — Parametric sweep (vary one geometric parameter):
%     results_table = trade_study('mode', 'sweep', ...
%                                 'sweep_param', 'Di', ...
%                                 'sweep_range', [0.015, 0.055], ...
%                                 'n_sweep', 8, ...)
%
%   Outputs:
%     results_table - MATLAB table with one row per configuration
%                     Columns: GrainType, Propellant, It_Ns, F_avg_N,
%                              F_max_N, Pc_max_MPa, tb_s, Isp_s, Class

    % ── Parse inputs ─────────────────────────────────────────────────
    ip = inputParser;
    addParameter(ip, 'mode',          'grain_compare');
    addParameter(ip, 'propellant',    'APCP_Aerotech');
    addParameter(ip, 'grain_type',    'BATES');
    addParameter(ip, 'Do',            0.075);
    addParameter(ip, 'Di',            0.025);
    addParameter(ip, 'L',             0.15);
    addParameter(ip, 'N',             4);
    addParameter(ip, 'At',            pi*(0.015/2)^2);
    addParameter(ip, 'epsilon',       8);
    addParameter(ip, 'P_amb',         101325);
    addParameter(ip, 'sweep_param',   'Di');
    addParameter(ip, 'sweep_range',   [0.015, 0.055]);
    addParameter(ip, 'n_sweep',       8);
    parse(ip, varargin{:});
    o = ip.Results;

    fprintf('\n====== Trade Study: %s ======\n\n', o.mode);

    all_results = {};

    switch o.mode

        % ─────────────────────────────────────────────────────────────
        case 'grain_compare'
        % Compare BATES, Star, WagonWheel with same propellant and nozzle
            prop   = propellant_db(o.propellant);
            nozzle = nozzle_model(prop, o.At, o.epsilon, o.P_amb);

            % BATES
            g1 = BATES(o.Do, o.Di, o.L, o.N, 500);
            r1 = ballistics_ode(g1, prop, nozzle);
            all_results{end+1} = r1;

            % Star grain (matched approximate web)
            Ri_star = o.Di/2 * 0.65;
            g2 = star_grain(o.Do, Ri_star, 0.45, 6, o.L * o.N, 500);
            r2 = ballistics_ode(g2, prop, nozzle);
            all_results{end+1} = r2;

            % Wagon wheel
            Rc    = o.Di/2 * 0.6;
            spk_w = o.Do * 0.06;
            spk_L = o.Do * 0.15;
            if Rc + spk_L < o.Do/2 * 0.9 && spk_w < 2*pi*Rc/6 * 0.8
                g3 = wagon_wheel(o.Do, Rc, 6, spk_w, spk_L, o.L*o.N, 500);
                r3 = ballistics_ode(g3, prop, nozzle);
                all_results{end+1} = r3;
            end

        % ─────────────────────────────────────────────────────────────
        case 'propellant_compare'
        % Compare all propellants on same grain
            propellants = {'APCP_Aerotech','APCP_HighPerf','KNSB','KNSU','KNDX'};
            g = BATES(o.Do, o.Di, o.L, o.N, 500);
            for i = 1:numel(propellants)
                try
                    prop   = propellant_db(propellants{i});
                    nozzle = nozzle_model(prop, o.At, o.epsilon, o.P_amb);
                    r      = ballistics_ode(g, prop, nozzle);
                    all_results{end+1} = r;
                catch ME
                    fprintf('Skipped %s: %s\n', propellants{i}, ME.message);
                end
            end

        % ─────────────────────────────────────────────────────────────
        case 'sweep'
        % Sweep one geometric parameter for BATES grain
            prop      = propellant_db(o.propellant);
            nozzle    = nozzle_model(prop, o.At, o.epsilon, o.P_amb);
            sweep_vec = linspace(o.sweep_range(1), o.sweep_range(2), o.n_sweep);

            for k = 1:numel(sweep_vec)
                Do = o.Do; Di = o.Di; L = o.L; N = o.N;
                switch o.sweep_param
                    case 'Di',  Di = sweep_vec(k);
                    case 'Do',  Do = sweep_vec(k);
                    case 'L',   L  = sweep_vec(k);
                    case 'N',   N  = round(sweep_vec(k));
                end
                if Di >= Do * 0.95, continue; end
                try
                    g = BATES(Do, Di, L, N, 400);
                    r = ballistics_ode(g, prop, nozzle);
                    r.sweep_val = sweep_vec(k);
                    all_results{end+1} = r;
                catch ME
                    fprintf('Sweep k=%d failed: %s\n', k, ME.message);
                end
            end

        otherwise
            error('Unknown mode: %s', o.mode);
    end

    % ── Build results table ───────────────────────────────────────────
    n = numel(all_results);
    GrainType   = cell(n,1);
    Propellant  = cell(n,1);
    It_Ns       = zeros(n,1);
    F_avg_N     = zeros(n,1);
    F_max_N     = zeros(n,1);
    Pc_max_MPa  = zeros(n,1);
    tb_s        = zeros(n,1);
    Isp_s       = zeros(n,1);
    Class       = cell(n,1);
    mp_kg       = zeros(n,1);

    for i = 1:n
        r = all_results{i};
        GrainType{i}   = r.grain.type;
        Propellant{i}  = r.prop.name;
        It_Ns(i)       = r.It;
        F_avg_N(i)     = r.F_avg;
        F_max_N(i)     = r.F_max;
        Pc_max_MPa(i)  = r.Pc_max / 1e6;
        tb_s(i)        = r.tb;
        Isp_s(i)       = r.Isp;
        Class{i}       = r.motor_class;
        mp_kg(i)       = r.mp;
    end

    results_table = table(GrainType, Propellant, It_Ns, F_avg_N, F_max_N, ...
                          Pc_max_MPa, tb_s, Isp_s, Class, mp_kg, ...
                          'VariableNames', ...
                          {'GrainType','Propellant','It_Ns','F_avg_N', ...
                           'F_max_N','Pc_max_MPa','tb_s','Isp_s','Class','mp_kg'});

    disp(results_table);

    % ── Comparison plots ──────────────────────────────────────────────
    plot_trade_study(all_results, o.mode, o.sweep_param);
end


function plot_trade_study(all_results, mode, sweep_param)
    n      = numel(all_results);
    colors = lines(n);

    figure('Name', sprintf('Trade Study — %s', mode), ...
           'Color', [1 1 1], 'Position', [80 80 1200 750]);

    % ── Panel 1: Thrust curves overlay ───────────────────────────────
    subplot(2,3,1);
    hold on; grid on; box on;
    for i = 1:n
        r   = all_results{i};
        lbl = sprintf('%s / %s', r.grain.type, r.prop.name);
        if strcmp(mode,'sweep')
            lbl = sprintf('%s=%.3f', sweep_param, r.sweep_val);
        end
        plot(r.t, r.F, 'Color', colors(i,:), 'LineWidth', 1.8, 'DisplayName', lbl);
    end
    xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust curves', 'FontWeight','normal');
    legend('Location','northeast','FontSize',7);

    % ── Panel 2: Pressure curves ──────────────────────────────────────
    subplot(2,3,2);
    hold on; grid on; box on;
    for i = 1:n
        r = all_results{i};
        plot(r.t, r.Pc/1e6, 'Color', colors(i,:), 'LineWidth', 1.8);
    end
    xlabel('Time [s]'); ylabel('Pc [MPa]');
    title('Chamber pressure', 'FontWeight','normal');

    % ── Panel 3: Total impulse bar chart ─────────────────────────────
    subplot(2,3,3);
    It_vals = cellfun(@(r) r.It, all_results);
    labels  = cellfun(@(r) sprintf('%s\n%s', r.grain.type, r.prop.name), ...
                      all_results, 'UniformOutput', false);
    if strcmp(mode,'sweep')
        labels = cellfun(@(r) sprintf('%.3f', r.sweep_val), all_results, 'UniformOutput', false);
    end
    bar(It_vals, 'FaceColor', 'flat', 'CData', colors(1:n,:));
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 20);
    ylabel('Total impulse [N·s]');
    title('Total impulse comparison', 'FontWeight','normal');
    grid on; box on;

    % ── Panel 4: Isp comparison ───────────────────────────────────────
    subplot(2,3,4);
    Isp_vals = cellfun(@(r) r.Isp, all_results);
    bar(Isp_vals, 'FaceColor', 'flat', 'CData', colors(1:n,:));
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 20);
    ylabel('Delivered I_{sp} [s]');
    title('Specific impulse comparison', 'FontWeight','normal');
    grid on; box on;

    % ── Panel 5: Peak pressure vs avg thrust scatter ──────────────────
    subplot(2,3,5);
    hold on; grid on; box on;
    for i = 1:n
        r = all_results{i};
        scatter(r.Pc_max/1e6, r.F_avg, 80, colors(i,:), 'filled');
        text(r.Pc_max/1e6, r.F_avg, sprintf('  %s', r.grain.type), 'FontSize',7);
    end
    xlabel('Peak chamber pressure [MPa]');
    ylabel('Average thrust [N]');
    title('Peak P_c vs. avg thrust', 'FontWeight','normal');

    % ── Panel 6: Burn duration ────────────────────────────────────────
    subplot(2,3,6);
    tb_vals = cellfun(@(r) r.tb, all_results);
    bar(tb_vals, 'FaceColor', 'flat', 'CData', colors(1:n,:));
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 20);
    ylabel('Burn duration [s]');
    title('Burn duration comparison', 'FontWeight','normal');
    grid on; box on;

    sgtitle(sprintf('Trade Study: %s', strrep(mode,'_',' ')), ...
            'FontSize', 13, 'FontWeight','normal');
end
