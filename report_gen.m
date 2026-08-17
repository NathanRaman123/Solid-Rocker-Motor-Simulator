%% report_gen.m — Auto-generate PDF summary report
%
%  Runs the full pipeline and saves all figures + a text summary
%  to a timestamped output folder. Use this to produce your final
%  "deliverable" — the kind of thing you'd show in a portfolio or
%  attach to a design review.
%
%  Usage:
%    report_gen()                     % uses default parameters
%    report_gen('propellant','KNSB')  % override propellant
%
%  Output: Creates folder  ./SRM_Report_YYYYMMDD_HHMMSS/
%          containing PNG plots + performance_summary.txt

function report_gen(varargin)

    ip = inputParser;
    addParameter(ip, 'propellant',  'APCP_Aerotech');
    addParameter(ip, 'Do',          0.075);
    addParameter(ip, 'Di',          0.025);
    addParameter(ip, 'L',           0.15);
    addParameter(ip, 'N',           4);
    addParameter(ip, 'At',          pi*(0.015/2)^2);
    addParameter(ip, 'epsilon',     8);
    parse(ip, varargin{:});
    o = ip.Results;

    % ── Create output folder ──────────────────────────────────────────
    timestamp  = datestr(now, 'yyyymmdd_HHMMSS');
    out_folder = sprintf('SRM_Report_%s', timestamp);
    mkdir(out_folder);
    fprintf('[Report] Output folder: %s\n', out_folder);

    % ── Run full simulation ───────────────────────────────────────────
    prop   = propellant_db(o.propellant);
    nozzle = nozzle_model(prop, o.At, o.epsilon, 101325);
    grain  = BATES(o.Do, o.Di, o.L, o.N, 500);
    res    = ballistics_ode(grain, prop, nozzle);

    % ── Figure 1: Grain cross-section ────────────────────────────────
    f1 = figure('Color',[1 1 1],'Visible','off','Position',[0 0 900 700]);
    grain_viz(grain);
    saveas(f1, fullfile(out_folder, '01_grain_geometry.png'));
    close(f1);

    % ── Figure 2: Thrust and pressure ────────────────────────────────
    f2 = figure('Color',[1 1 1],'Visible','off','Position',[0 0 1000 500]);
    subplot(1,2,1);
    plot(res.t, res.F, 'b-', 'LineWidth', 2.5);
    hold on; grid on; box on;
    xline(res.tb, 'r--', 'Burnout', 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('Thrust [N]');
    title(sprintf('Thrust — Class %s  (%.1f N·s)', res.motor_class, res.It), ...
          'FontWeight','normal');

    subplot(1,2,2);
    plot(res.t, res.Pc/1e6, 'r-', 'LineWidth', 2.5);
    hold on; grid on; box on;
    xline(res.tb, 'r--', 'Burnout', 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel('Pressure [MPa]');
    title(sprintf('Chamber pressure (peak %.2f MPa)', res.Pc_max/1e6), ...
          'FontWeight','normal');
    sgtitle(sprintf('%s / %s — Performance', grain.type, prop.name), ...
            'FontSize',12,'FontWeight','normal');
    saveas(f2, fullfile(out_folder, '02_thrust_pressure.png'));
    close(f2);

    % ── Figure 3: Trade study — grain comparison ──────────────────────
    T = trade_study('mode', 'grain_compare', ...
                    'propellant', o.propellant, ...
                    'Do', o.Do, 'Di', o.Di, 'L', o.L, 'N', o.N, ...
                    'At', o.At, 'epsilon', o.epsilon);
    f3 = gcf;
    saveas(f3, fullfile(out_folder, '03_grain_trade_study.png'));
    close(f3);

    % ── Figure 4: Propellant comparison ──────────────────────────────
    T2 = trade_study('mode', 'propellant_compare', ...
                     'Do', o.Do, 'Di', o.Di, 'L', o.L, 'N', o.N, ...
                     'At', o.At, 'epsilon', o.epsilon);
    f4 = gcf;
    saveas(f4, fullfile(out_folder, '04_propellant_trade_study.png'));
    close(f4);

    % ── Figure 5: Validation ─────────────────────────────────────────
    val = validation(res, 'builtin');
    f5  = gcf;
    saveas(f5, fullfile(out_folder, '05_validation.png'));
    close(f5);

    % ── Write text summary ────────────────────────────────────────────
    fid = fopen(fullfile(out_folder, 'performance_summary.txt'), 'w');
    fprintf(fid, 'Solid Rocket Motor Simulator — Performance Summary\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, '%s\n\n', repmat('=', 1, 52));

    fprintf(fid, 'DESIGN INPUTS\n');
    fprintf(fid, '  Grain type     : %s\n', grain.type);
    fprintf(fid, '  OD             : %.1f mm\n', o.Do*1e3);
    fprintf(fid, '  Core dia       : %.1f mm\n', o.Di*1e3);
    fprintf(fid, '  Segment length : %.1f mm\n', o.L*1e3);
    fprintf(fid, '  Segments       : %d\n', o.N);
    fprintf(fid, '  Propellant     : %s\n', prop.name);
    fprintf(fid, '  Throat dia     : %.1f mm\n', nozzle.Dt*1e3);
    fprintf(fid, '  Expansion ratio: %.1f\n\n', o.epsilon);

    fprintf(fid, 'PERFORMANCE RESULTS\n');
    fprintf(fid, '  Motor class    : %s\n', res.motor_class);
    fprintf(fid, '  Total impulse  : %.1f N·s\n', res.It);
    fprintf(fid, '  Burn duration  : %.3f s\n', res.tb);
    fprintf(fid, '  Peak thrust    : %.1f N\n', res.F_max);
    fprintf(fid, '  Avg thrust     : %.1f N\n', res.F_avg);
    fprintf(fid, '  Peak Pc        : %.2f MPa\n', res.Pc_max/1e6);
    fprintf(fid, '  Delivered Isp  : %.1f s\n', res.Isp);
    fprintf(fid, '  Propellant mass: %.3f kg\n\n', res.mp);

    fprintf(fid, 'VALIDATION (vs J350W builtin reference)\n');
    fprintf(fid, '  R²             : %.4f\n', val.r2);
    fprintf(fid, '  Total It error : %.1f%%\n', val.err_It);
    fprintf(fid, '  Peak F error   : %.1f%%\n', val.err_Fpeak);
    fprintf(fid, '  RMSE           : %.2f N\n', val.rmse);
    fclose(fid);

    fprintf('\n[Report] Complete. Files saved to: ./%s/\n', out_folder);
    fprintf('[Report] Contents:\n');
    fprintf('   01_grain_geometry.png\n');
    fprintf('   02_thrust_pressure.png\n');
    fprintf('   03_grain_trade_study.png\n');
    fprintf('   04_propellant_trade_study.png\n');
    fprintf('   05_validation.png\n');
    fprintf('   performance_summary.txt\n\n');
end
