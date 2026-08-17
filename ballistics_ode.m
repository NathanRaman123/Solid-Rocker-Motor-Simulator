function results = ballistics_ode(grain, prop, nozzle, varargin)
% BALLISTICS_ODE  Internal ballistics ODE solver for solid rocket motors
%
%   Integrates the chamber pressure ODE using MATLAB's ode45 (RK4/5).
%   The governing equation is a mass balance on the chamber gas:
%
%     dPc/dt = (rho_prop * Ab(w) * r(Pc) * cstar - Pc * At) * cstar / Vc(w)
%
%   where:
%     Pc   = chamber pressure [Pa]
%     w    = web regression [m]  (integrated from dw/dt = r)
%     Ab   = burn area as function of web [m^2]
%     r    = burn rate = a * Pc^n  [m/s]
%     At   = throat area [m^2]
%     Vc   = chamber free volume [m^3]
%     cstar= characteristic velocity [m/s]
%     rho  = propellant density [kg/m^3]
%
%   State vector: y = [Pc; w]
%
%   Inputs:
%     grain    - Grain geometry struct from BATES/star_grain/wagon_wheel
%     prop     - Propellant struct from propellant_db
%     nozzle   - Nozzle struct from nozzle_model
%
%   Optional name-value pairs:
%     'Pc_init'   - Initial chamber pressure [Pa]  (default: 101325 Pa = 1 atm)
%     'V_casing'  - Total casing internal volume [m^3]  (default: computed from grain)
%     'dt_max'    - Max ODE timestep [s]  (default: 1e-4)
%     'plot'      - true/false: plot results immediately  (default: false)
%
%   Outputs:
%     results struct:
%       .t          time vector [s]
%       .Pc         chamber pressure [Pa]
%       .F          thrust [N]
%       .w          web regression [m]
%       .Ab         burn area [m^2]
%       .r          burn rate [m/s]
%       .It         total impulse [N·s]
%       .Isp        delivered specific impulse [s]
%       .tb         burn duration [s]
%       .Pc_max     peak chamber pressure [Pa]
%       .F_avg      average thrust [N]
%       .F_max      peak thrust [N]
%       .motor_class NAR/TRA letter classification
%       .mp         propellant mass [kg]
%
%   Example:
%     g = BATES(0.075, 0.025, 0.15, 4, 500);
%     p = propellant_db('APCP_Aerotech');
%     n = nozzle_model(p, pi*(0.015/2)^2, 8, 101325);
%     res = ballistics_ode(g, p, n, 'plot', true);

    % ── Input validation ──────────────────────────────────────────────
    narginchk(3, Inf);

    % ── Parse options ─────────────────────────────────────────────────
    ip = inputParser;
    addParameter(ip, 'Pc_init',  101325);
    addParameter(ip, 'V_casing', []);
    addParameter(ip, 'dt_max',   1e-4);
    addParameter(ip, 'plot',     false);
    parse(ip, varargin{:});

    Pc0      = ip.Results.Pc_init;
    dt_max   = ip.Results.dt_max;
    do_plot  = ip.Results.plot;

    % ── Interpolation functions for Ab and Vp ────────────────────────
    % The ODE needs Ab(w) and Vp(w) — create interpolants from grain data
    Ab_interp = @(w) interp1(grain.w, grain.Ab, ...
                             min(max(w, 0), grain.w_max), 'linear', 0);
    Vp_interp = @(w) interp1(grain.w, grain.Vp, ...
                             min(max(w, 0), grain.w_max), 'linear', 0);

    % ── Chamber free volume ───────────────────────────────────────────
    % V_chamber = V_casing_total - V_propellant_remaining
    % Initial propellant volume (t=0)
    Vp0 = Vp_interp(0);
    % Estimate casing volume: propellant volume + 15% ullage
    if isempty(ip.Results.V_casing)
        V_casing = Vp0 * 1.15;
    else
        V_casing = ip.Results.V_casing;
    end
    % Free volume as function of web
    Vc_func = @(w) V_casing - Vp_interp(w);

    % ── ODE definition ────────────────────────────────────────────────
    % State: y(1) = Pc [Pa],  y(2) = w [m]
    function dydt = odes(~, y)
        Pc_now = max(y(1), 101325);   % clamp to ambient
        w_now  = min(max(y(2), 0), grain.w_max);

        % Burn rate and burn area
        r_now  = prop.a * Pc_now^prop.n;
        Ab_now = Ab_interp(w_now);
        Vc_now = Vc_func(w_now);

        % After burnout, no more mass generation
        if w_now >= grain.w_max || Ab_now <= 0
            r_now  = 0;
            Ab_now = 0;
        end

        % Mass generation rate: mdot_gen = rho * Ab * r
        % Mass loss through nozzle: mdot_out = Pc * At / cstar
        % dPc/dt = cstar/Vc * (rho*Ab*r - Pc*At/cstar) * (R_gas*Tc/(gamma*cstar))
        % Simplified quasi-steady form (ref: Sutton RPE ch.12):
        dPc_dt = (prop.cstar / Vc_now) * ...
                 (prop.rho * Ab_now * r_now - Pc_now * nozzle.At / prop.cstar);

        % Web regression rate
        dw_dt = r_now;

        dydt = [dPc_dt; dw_dt];
    end

    % ── Event: burnout ────────────────────────────────────────────────
    function [val, isterminal, direction] = burnout_event(~, y)
        val        = grain.w_max - y(2) - 1e-6;  % zero when w → w_max
        isterminal = 0;   % don't stop — let pressure decay naturally
        direction  = -1;
    end

    function [val, isterminal, direction] = pressure_floor_event(~, y)
        val        = y(1) - 1.05 * 101325;  % pressure back to ~ambient
        isterminal = 1;                       % stop integration
        direction  = -1;
    end

    % ── ODE solve ─────────────────────────────────────────────────────
    y0       = [Pc0; 0];
    t_span   = [0, 60];    % generous upper bound — event will stop early

    opts = odeset('Events',   @pressure_floor_event, ...
                  'MaxStep',  dt_max, ...
                  'RelTol',   1e-6, ...
                  'AbsTol',   1e-8);

    fprintf('[ODE] Integrating ballistics...\n');
    [t_raw, y_raw] = ode45(@odes, t_span, y0, opts);

    % ── Post-process ──────────────────────────────────────────────────
    Pc_raw = y_raw(:,1);
    w_raw  = y_raw(:,2);

    % Derived quantities at each time step
    Ab_raw = arrayfun(@(w) Ab_interp(w), w_raw);
    r_raw  = prop.a .* Pc_raw.^prop.n;
    r_raw(w_raw >= grain.w_max) = 0;
    Ab_raw(w_raw >= grain.w_max) = 0;

    % Thrust: F = Cf * At * Pc
    Cf_raw = arrayfun(@(Pc) nozzle.Cf(Pc), Pc_raw);
    F_raw  = Cf_raw .* nozzle.At .* Pc_raw;
    F_raw  = max(F_raw, 0);

    % Burn duration (when web = w_max)
    idx_burnout = find(w_raw >= grain.w_max, 1, 'first');
    if isempty(idx_burnout)
        idx_burnout = length(t_raw);
    end
    tb = t_raw(idx_burnout);

    % Total impulse (trapz integration under F vs t)
    It = trapz(t_raw, F_raw);

    % Propellant mass
    mp = prop.rho * Vp_interp(0);

    % Delivered Isp
    Isp_del = It / (mp * 9.80665);

    % Motor classification
    motor_class = classify_motor(It);

    % ── Package results ───────────────────────────────────────────────
    results.t           = t_raw;
    results.Pc          = Pc_raw;
    results.F           = F_raw;
    results.Cf          = Cf_raw;
    results.w           = w_raw;
    results.Ab          = Ab_raw;
    results.r           = r_raw;
    results.It          = It;
    results.Isp         = Isp_del;
    results.tb          = tb;
    results.Pc_max      = max(Pc_raw);
    results.F_avg       = mean(F_raw(F_raw > 0.01*max(F_raw)));
    results.F_max       = max(F_raw);
    results.motor_class = motor_class;
    results.mp          = mp;
    results.grain       = grain;
    results.prop        = prop;
    results.nozzle      = nozzle;

    % ── Console summary ───────────────────────────────────────────────
    fprintf('\n========== Motor Performance Summary ==========\n');
    fprintf('  Grain type    : %s (%s)\n', grain.type, grain.profile_type);
    fprintf('  Propellant    : %s\n', prop.name);
    fprintf('  Total impulse : %.1f N·s  → Class %s\n', It, motor_class);
    fprintf('  Burn duration : %.3f s\n', tb);
    fprintf('  Peak pressure : %.2f MPa\n', max(Pc_raw)/1e6);
    fprintf('  Avg pressure  : %.2f MPa\n', mean(Pc_raw(1:idx_burnout))/1e6);
    fprintf('  Peak thrust   : %.1f N\n', max(F_raw));
    fprintf('  Avg thrust    : %.1f N\n', results.F_avg);
    fprintf('  Delivered Isp : %.1f s\n', Isp_del);
    fprintf('  Prop mass     : %.3f kg\n', mp);
    fprintf('===============================================\n\n');

    % ── Optional plot ─────────────────────────────────────────────────
    if do_plot
        plot_ballistics(results);
    end
end


% ── Motor classification (NAR/TRA standard) ───────────────────────────
function cls = classify_motor(It)
    % Each letter class spans 1.26–2.50 N·s for A, doubling each step
    boundaries = [0, 1.25, 2.5, 5, 10, 20, 40, 80, 160, 320, 640, ...
                  1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840];
    letters    = 'ABCDEFGHIJKLMNOPQRS';
    idx = find(It > boundaries, 1, 'last');
    if isempty(idx) || idx > length(letters)
        cls = '?';
    else
        cls = letters(idx);
    end
end


% ── Internal plotting ─────────────────────────────────────────────────
function plot_ballistics(res)
    figure('Name', sprintf('Ballistics — %s / %s', res.grain.type, res.prop.name), ...
           'Color', [1 1 1], 'Position', [100 100 1100 700]);

    t  = res.t;
    tb = res.tb;

    % Panel 1: Thrust
    subplot(2,3,1);
    plot(t, res.F, 'b-', 'LineWidth', 2);
    hold on; grid on; box on;
    xline(tb, 'r--', 'LineWidth', 1, 'DisplayName', 'Burnout');
    xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust F(t)', 'FontWeight', 'normal');
    legend(sprintf('It=%.1fNs  Class %s', res.It, res.motor_class), ...
           'Location', 'northeast');

    % Panel 2: Chamber pressure
    subplot(2,3,2);
    plot(t, res.Pc/1e6, 'r-', 'LineWidth', 2);
    hold on; grid on; box on;
    xline(tb, 'r--', 'LineWidth', 1);
    yline(res.Pc_max/1e6, 'k:', 'LineWidth', 1);
    xlabel('Time [s]'); ylabel('Chamber pressure [MPa]');
    title('Chamber Pressure P_c(t)', 'FontWeight', 'normal');
    text(0.98, 0.95, sprintf('Peak: %.2f MPa', res.Pc_max/1e6), ...
         'Units','normalized','HorizontalAlignment','right','FontSize',8);

    % Panel 3: Burn area
    subplot(2,3,3);
    plot(res.w*1e3, res.Ab*1e4, 'm-', 'LineWidth', 2);
    grid on; box on;
    xlabel('Web regression [mm]'); ylabel('Burn area [cm²]');
    title('Burn Area A_b(w)', 'FontWeight', 'normal');

    % Panel 4: Burn rate vs pressure
    subplot(2,3,4);
    P_range = linspace(0.1e6, res.Pc_max*1.1, 300);
    r_range = res.prop.a .* P_range.^res.prop.n;
    plot(P_range/1e6, r_range*1e3, 'k-', 'LineWidth', 2);
    hold on; grid on; box on;
    scatter(res.Pc/1e6, res.r*1e3, 4, t, 'filled');
    colorbar; xlabel('Pressure [MPa]'); ylabel('Burn rate [mm/s]');
    title(sprintf('Saint-Robert''s Law (n=%.3f)', res.prop.n), 'FontWeight','normal');

    % Panel 5: Web vs time
    subplot(2,3,5);
    plot(t, res.w*1e3, 'g-', 'LineWidth', 2);
    hold on; grid on; box on;
    yline(res.grain.w_max*1e3, 'k--', 'w_{max}', 'LineWidth', 1);
    xlabel('Time [s]'); ylabel('Web regression [mm]');
    title('Web Regression w(t)', 'FontWeight', 'normal');

    % Panel 6: Performance summary text
    subplot(2,3,6);
    axis off;
    summary = {
        sprintf('Motor Class:    %s', res.motor_class);
        sprintf('Total Impulse:  %.1f N·s', res.It);
        sprintf('Burn Duration:  %.3f s', res.tb);
        sprintf('Peak Thrust:    %.1f N', res.F_max);
        sprintf('Avg Thrust:     %.1f N', res.F_avg);
        sprintf('Peak Pc:        %.2f MPa', res.Pc_max/1e6);
        sprintf('Delivered Isp:  %.1f s', res.Isp);
        sprintf('Prop Mass:      %.3f kg', res.mp);
        '';
        sprintf('Grain:          %s', res.grain.type);
        sprintf('Propellant:     %s', res.prop.name);
        sprintf('Nozzle Dt:      %.1f mm', res.nozzle.Dt*1e3);
        sprintf('Expansion ε:    %.1f', res.nozzle.epsilon);
    };
    for i = 1:numel(summary)
        text(0.05, 1.0 - (i-1)*0.075, summary{i}, ...
             'Units','normalized','FontSize',9,'FontName','Monospaced');
    end
    title('Performance Summary', 'FontWeight','normal');

    sgtitle(sprintf('Internal Ballistics — %s / %s', ...
            res.grain.type, res.prop.name), ...
            'FontSize', 13, 'FontWeight', 'normal');
end
