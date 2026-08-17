function [r, r_curve] = burn_rate(prop, P, varargin)
% BURN_RATE  Saint-Robert's (Vielle's) burn rate law
%
%   Computes instantaneous burn rate using:
%       r = a * P^n    [m/s]
%
%   where P is chamber pressure [Pa], a and n are propellant-specific
%   coefficients from propellant_db().
%
%   Inputs:
%     prop        - Propellant struct from propellant_db()
%     P           - Chamber pressure, scalar or vector [Pa]
%
%   Optional name-value pairs:
%     'plot'      - true/false: plot r vs P curve (default false)
%     'P_range'   - [P_min P_max] for plot [Pa]  (default [0.5 10] MPa)
%     'T_ambient' - Ambient temperature [K] for temperature sensitivity
%                   (default 294 K — standard day)
%
%   Outputs:
%     r           - Burn rate at input pressure(s) [m/s]
%     r_curve     - Struct with full pressure sweep:
%                     .P      pressure vector [Pa]
%                     .r      burn rate vector [m/s]
%                     .sigma  temperature sensitivity [1/K]  (estimated)
%
%   Stability check:
%     A motor is stable when n < 1. If n >= 1, burn rate increases
%     faster than the nozzle can exhaust mass — runaway pressurization.
%     This function warns if n >= 0.8 (approaching instability).
%
%   Example:
%     prop = propellant_db('APCP_Aerotech');
%     r = burn_rate(prop, 3.5e6);          % at 3.5 MPa
%     [~, rc] = burn_rate(prop, [], 'plot', true);

    p = inputParser;
    addOptional(p, 'plot_flag', false);
    addParameter(p, 'plot', false);
    addParameter(p, 'P_range', [0.5e6, 12e6]);
    addParameter(p, 'T_ambient', 294);
    parse(p, varargin{:});
    do_plot   = p.Results.plot || p.Results.plot_flag;
    P_range   = p.Results.P_range;
    T_amb     = p.Results.T_ambient;

    % ── Stability warning ─────────────────────────────────────────────
    if prop.n >= 1.0
        warning('BURN_RATE:unstable', ...
            'Pressure exponent n=%.3f >= 1.0 — UNSTABLE motor! Runaway pressurization risk.', prop.n);
    elseif prop.n >= 0.8
        warning('BURN_RATE:nearUnstable', ...
            'Pressure exponent n=%.3f >= 0.8 — approaching instability.', prop.n);
    end

    % ── Burn rate at requested pressure(s) ───────────────────────────
    if isempty(P)
        r = [];
    else
        r = prop.a .* P.^prop.n;
    end

    % ── Full pressure sweep for curve ────────────────────────────────
    P_sweep = linspace(P_range(1), P_range(2), 500);
    r_sweep = prop.a .* P_sweep.^prop.n;

    % Temperature sensitivity coefficient (sigma)
    % Approximate: sigma ≈ (1/r) * dr/dT ≈ n * (da/dT)/a
    % Using empirical estimate: ~0.2%/K for APCP, ~0.3%/K for KN propellants
    if contains(prop.name, 'KN', 'IgnoreCase', true)
        sigma_r = 0.003;   % 1/K
    else
        sigma_r = 0.002;
    end

    % Temperature-corrected burn rate (if not standard day)
    T_std   = 294;   % K
    dT      = T_amb - T_std;
    r_T_corrected = r_sweep .* exp(sigma_r * dT);

    r_curve.P               = P_sweep;
    r_curve.r               = r_sweep;
    r_curve.r_cold          = r_sweep .* exp(sigma_r * (-30));  % -30 K cold day
    r_curve.r_hot           = r_sweep .* exp(sigma_r * (+30));  % +30 K hot day
    r_curve.r_T_corrected   = r_T_corrected;
    r_curve.sigma           = sigma_r;
    r_curve.T_amb           = T_amb;

    % ── Optional plot ─────────────────────────────────────────────────
    if do_plot
        figure('Name', sprintf('Burn Rate — %s', prop.name), ...
               'Color', [1 1 1], 'Position', [150 150 750 450]);
        hold on; grid on; box on;

        % Temperature sensitivity band
        fill([P_sweep, fliplr(P_sweep)] / 1e6, ...
             [r_curve.r_cold, fliplr(r_curve.r_hot)] * 1e3, ...
             [0.8 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
             'DisplayName', '±30K temperature band');

        % Nominal curve
        plot(P_sweep/1e6, r_sweep*1e3, 'b-', 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%s (T=%.0fK)', prop.name, T_std));

        % Mark if specific pressure was requested
        if ~isempty(P) && isscalar(P)
            plot(P/1e6, r*1e3, 'ro', 'MarkerSize', 8, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Query: P=%.2f MPa → r=%.2f mm/s', P/1e6, r*1e3));
        end

        % Slope annotation (n value shown on log-log is slope)
        xlabel('Chamber pressure [MPa]');
        ylabel('Burn rate [mm/s]');
        title(sprintf('Saint-Robert''s Law: r = a·P^n  (a=%.2e, n=%.3f)', prop.a, prop.n), ...
              'FontWeight', 'normal');
        legend('Location', 'northwest');

        % Inset: log-log to show linearity
        ax2 = axes('Position', [0.60 0.18 0.28 0.30]);
        loglog(ax2, P_sweep/1e6, r_sweep*1e3, 'b-', 'LineWidth', 1.5);
        grid(ax2, 'on');
        xlabel(ax2, 'P [MPa]', 'FontSize', 8);
        ylabel(ax2, 'r [mm/s]', 'FontSize', 8);
        title(ax2, 'Log-log (slope = n)', 'FontSize', 8, 'FontWeight', 'normal');
    end

    if ~isempty(P) && isscalar(P)
        fprintf('[BurnRate] %s | P=%.2f MPa → r=%.3f mm/s\n', ...
            prop.name, P/1e6, r*1e3);
    end
end
