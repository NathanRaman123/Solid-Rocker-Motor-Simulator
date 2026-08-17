function val = validation(results, eng_file_path, varargin)
% VALIDATION  Compare simulator output against published .eng thrust data
%
%   Loads a RASP/OpenRocket .eng file (standard format used by AeroTech,
%   Cesaroni, and the amateur rocketry community) and computes error
%   metrics between your simulation and the certified data.
%
%   .eng file format (text):
%     Line 1: ; comment
%     Line 2: <MotorName> <Dia_mm> <Len_mm> <delays> <mp_g> <total_g> <mfr>
%     Lines 3+: <time_s> <thrust_N>   (space-separated pairs)
%     Last line: 0.000 0.000
%
%   Inputs:
%     results       - Output struct from ballistics_ode()
%     eng_file_path - Path to .eng file, OR 'builtin' to use a
%                     hardcoded reference curve for testing
%
%   Optional:
%     'motor_name'  - String label for plot title
%     'plot'        - true/false (default true)
%
%   Outputs:
%     val struct:
%       .t_ref      reference time [s]
%       .F_ref      reference thrust [N]
%       .t_sim      simulation time [s]
%       .F_sim      simulation thrust (resampled to ref) [N]
%       .It_ref     reference total impulse [N·s]
%       .It_sim     simulation total impulse [N·s]
%       .err_It     total impulse error [%]
%       .err_Fpeak  peak thrust error [%]
%       .rmse       RMSE of thrust curve [N]
%       .r2         R² correlation coefficient
%
%   Example:
%     res = ballistics_ode(grain, prop, nozzle);
%     val = validation(res, 'AeroTech_J350W.eng', 'plot', true);
%
%   Finding .eng files:
%     http://www.thrustcurve.org  (search any certified motor, download .eng)
%     https://www.rocketmotorpro.com

    ip = inputParser;
    addParameter(ip, 'motor_name', 'Reference motor');
    addParameter(ip, 'plot', true);
    parse(ip, varargin{:});
    motor_name = ip.Results.motor_name;
    do_plot    = ip.Results.plot;

    % ── Load reference data ───────────────────────────────────────────
    if strcmp(eng_file_path, 'builtin')
        % Built-in reference: synthetic J-class thrust curve
        % Approximates AeroTech J350W profile for testing
        [t_ref, F_ref, meta] = builtin_reference_curve();
        motor_name = 'AeroTech J350W (builtin reference)';
    else
        [t_ref, F_ref, meta] = parse_eng_file(eng_file_path);
        if isempty(motor_name)
            motor_name = meta.name;
        end
    end

    fprintf('[Validation] Reference motor: %s\n', motor_name);
    fprintf('[Validation] Reference: It=%.1fNs  Fpk=%.0fN  tb=%.2fs\n', ...
        trapz(t_ref, F_ref), max(F_ref), t_ref(end));

    % ── Resample simulation to reference time grid ────────────────────
    t_sim_raw = results.t;
    F_sim_raw = results.F;

    % Extend sim to reference duration if needed (pad with zeros)
    t_end = max(t_ref(end), t_sim_raw(end));
    if t_sim_raw(end) < t_ref(end)
        t_sim_raw = [t_sim_raw; t_ref(end)];
        F_sim_raw = [F_sim_raw; 0];
    end

    F_sim_interp = interp1(t_sim_raw, F_sim_raw, t_ref, 'linear', 0);
    F_sim_interp = max(F_sim_interp, 0);

    % ── Error metrics ─────────────────────────────────────────────────
    It_ref   = trapz(t_ref, F_ref);
    It_sim   = trapz(t_ref, F_sim_interp);
    err_It   = (It_sim - It_ref) / It_ref * 100;

    Fpeak_ref = max(F_ref);
    Fpeak_sim = max(F_sim_interp);
    err_Fpeak = (Fpeak_sim - Fpeak_ref) / Fpeak_ref * 100;

    % RMSE
    rmse = sqrt(mean((F_sim_interp - F_ref).^2));

    % R² coefficient
    SS_res = sum((F_ref - F_sim_interp).^2);
    SS_tot = sum((F_ref - mean(F_ref)).^2);
    r2     = 1 - SS_res / SS_tot;

    % ── Console output ────────────────────────────────────────────────
    fprintf('\n--- Validation Results ---\n');
    fprintf('  Total impulse: Ref=%.1f Ns  Sim=%.1f Ns  Error=%.1f%%\n', ...
        It_ref, It_sim, err_It);
    fprintf('  Peak thrust:   Ref=%.0f N   Sim=%.0f N   Error=%.1f%%\n', ...
        Fpeak_ref, Fpeak_sim, err_Fpeak);
    fprintf('  RMSE:          %.2f N\n', rmse);
    fprintf('  R²:            %.4f\n', r2);
    fprintf('--------------------------\n\n');

    % ── Package output ────────────────────────────────────────────────
    val.t_ref     = t_ref;
    val.F_ref     = F_ref;
    val.t_sim     = t_ref;
    val.F_sim     = F_sim_interp;
    val.It_ref    = It_ref;
    val.It_sim    = It_sim;
    val.err_It    = err_It;
    val.err_Fpeak = err_Fpeak;
    val.rmse      = rmse;
    val.r2        = r2;
    val.meta      = meta;
    val.motor_name = motor_name;

    % ── Plot ──────────────────────────────────────────────────────────
    if do_plot
        figure('Name', sprintf('Validation — %s', motor_name), ...
               'Color', [1 1 1], 'Position', [120 120 950 550]);

        subplot(1,2,1);
        hold on; grid on; box on;
        plot(t_ref, F_ref, 'k-', 'LineWidth', 2.5, 'DisplayName', ...
             sprintf('Reference: %s', motor_name));
        plot(t_ref, F_sim_interp, 'b--', 'LineWidth', 2, 'DisplayName', ...
             sprintf('Simulation (R²=%.3f)', r2));
        % Error shading
        fill([t_ref; flipud(t_ref)], [F_ref; flipud(F_sim_interp)], ...
             [1 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4,'HandleVisibility','off');
        xlabel('Time [s]'); ylabel('Thrust [N]');
        title('Thrust curve: simulation vs. reference', 'FontWeight','normal');
        legend('Location','northeast');

        % Annotate metrics
        ax = gca;
        text(0.97, 0.50, sprintf('It err: %.1f%%\nFpk err: %.1f%%\nRMSE: %.1f N\nR²: %.4f', ...
             err_It, err_Fpeak, rmse, r2), ...
             'Units','normalized','HorizontalAlignment','right', ...
             'FontSize',9,'BackgroundColor',[0.95 0.95 0.95], ...
             'EdgeColor',[0.8 0.8 0.8]);

        subplot(1,2,2);
        residuals = F_sim_interp - F_ref;
        plot(t_ref, residuals, 'r-', 'LineWidth', 1.5);
        hold on; grid on; box on;
        yline(0, 'k-', 'LineWidth', 1);
        yline(rmse, 'b--', 'RMSE', 'LineWidth', 1);
        yline(-rmse, 'b--', 'LineWidth', 1);
        xlabel('Time [s]'); ylabel('Residual F_{sim} - F_{ref} [N]');
        title('Thrust residuals', 'FontWeight','normal');

        sgtitle(sprintf('Validation: %s', motor_name), ...
                'FontSize',12,'FontWeight','normal');
    end
end


% ── Parse a standard RASP .eng file ──────────────────────────────────
function [t, F, meta] = parse_eng_file(filepath)
    fid = fopen(filepath, 'r');
    if fid < 0
        error('Cannot open .eng file: %s', filepath);
    end

    meta = struct('name','','diameter_mm',0,'length_mm',0,'mp_g',0);
    t = []; F = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || line(1) == ';'
            continue;   % Skip comments
        end

        parts = strsplit(line);
        if numel(parts) >= 6 && ~all(isstrprop(parts{1}(1), 'digit'))
            % Header line: Name Dia Len Delays Mp_g Mtot_g Mfr
            meta.name        = parts{1};
            meta.diameter_mm = str2double(parts{2});
            meta.length_mm   = str2double(parts{3});
            meta.mp_g        = str2double(parts{5});
        elseif numel(parts) >= 2 && all(isstrprop(parts{1}(1), 'digit') | parts{1}(1) == '.')
            % Data line: time thrust
            t_val = str2double(parts{1});
            F_val = str2double(parts{2});
            if ~isnan(t_val) && ~isnan(F_val)
                t(end+1) = t_val; %#ok<AGROW>
                F(end+1) = F_val; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    t = t(:); F = F(:);

    % Ensure starts at t=0
    if t(1) > 0
        t = [0; t]; F = [0; F];
    end
    % Ensure ends at 0 thrust
    if F(end) > 0
        t = [t; t(end)+0.01]; F = [F; 0];
    end
end


% ── Built-in reference curve (J350W approximation) ───────────────────
function [t, F, meta] = builtin_reference_curve()
    % Approximate AeroTech J350W certified data
    % Based on publicly available thrustcurve.org data
    t_data = [0.000, 0.010, 0.030, 0.060, 0.100, 0.160, 0.230, 0.320, ...
              0.430, 0.560, 0.700, 0.850, 1.000, 1.150, 1.280, 1.380, ...
              1.440, 1.480, 1.510, 1.530, 1.540];
    F_data = [0,    280,   370,   390,   400,   395,   385,   375, ...
              368,  360,   355,   350,   345,   340,   330,   310, ...
              260,  180,    80,    20,     0];

    % Smooth slightly
    t = linspace(0, t_data(end), 200)';
    F = max(interp1(t_data, F_data, t, 'pchip', 0), 0);

    meta.name        = 'J350W';
    meta.diameter_mm = 54;
    meta.length_mm   = 325;
    meta.mp_g        = 242;
end
