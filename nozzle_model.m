function nozzle = nozzle_model(prop, At, expansion_ratio, P_amb)
% NOZZLE_MODEL  Isentropic nozzle flow for solid rocket motors
%
%   Inputs:
%     prop            - Propellant struct from propellant_db()
%     At              - Throat area [m^2]
%     expansion_ratio - Ae/At [-]
%     P_amb           - Ambient pressure [Pa]
%
%   Outputs:
%     nozzle struct including .Cf function handle (Pc -> Cf scalar)

    gamma = prop.gamma;
    R_u   = 8.314;
    R_gas = R_u / prop.M;
    Tc    = prop.Tc;

    % ── Geometry ──────────────────────────────────────────────────────
    Ae = At * expansion_ratio;
    Dt = 2 * sqrt(At / pi);
    De = 2 * sqrt(Ae / pi);

    % ── c* ────────────────────────────────────────────────────────────
    cstar_theory = sqrt(gamma * R_gas * Tc) / ...
                   (gamma * sqrt((2/(gamma+1))^((gamma+1)/(gamma-1))));
    cstar = prop.cstar;

    % ── Exit Mach — solve ONCE here, never again ──────────────────────
    Me = solve_exit_mach(gamma, expansion_ratio);

    % ── Design-point exit pressure and velocity ───────────────────────
    Pc_design = 3.5e6;
    Pe_design = Pc_design / (1 + (gamma-1)/2 * Me^2)^(gamma/(gamma-1));
    Ve_ideal  = Me * sqrt(gamma * R_gas * Tc / (1 + (gamma-1)/2 * Me^2));

    % ── Cf function — Me is captured in closure, NO fsolve inside ─────
    Cf_func = @(Pc) cf_from_me(Pc, P_amb, gamma, Me, expansion_ratio);
    nozzle.Cf = Cf_func;

    % ── Package ───────────────────────────────────────────────────────
    nozzle.At           = At;
    nozzle.Ae           = Ae;
    nozzle.Dt           = Dt;
    nozzle.De           = De;
    nozzle.epsilon      = expansion_ratio;
    nozzle.Me           = Me;
    nozzle.cstar        = cstar;
    nozzle.cstar_theory = cstar_theory;
    nozzle.Ve_ideal     = Ve_ideal;
    nozzle.Pe_design    = Pe_design;
    nozzle.Isp_vac      = prop.Isp_vac;
    nozzle.gamma        = gamma;
    nozzle.R_gas        = R_gas;
    nozzle.Tc           = Tc;
    nozzle.P_amb        = P_amb;
    nozzle.Pc_design    = Pc_design;

    fprintf('[Nozzle] Dt=%.1fmm  De=%.1fmm  e=%.1f  Me=%.2f  c*=%.0fm/s  Isp(vac)=%ds\n', ...
        Dt*1e3, De*1e3, expansion_ratio, Me, cstar, prop.Isp_vac);
    fprintf('         Pe(design@%.1fMPa)=%.1fkPa  Ve=%.0fm/s\n', ...
        Pc_design/1e6, Pe_design/1e3, Ve_ideal);
end


% ── Solve exit Mach from area ratio (bisection — no toolbox needed) ───
function Me = solve_exit_mach(gamma, epsilon)
    % Supersonic root of: A/At = f(M)
    % Bracket [1.001, 15] always contains exactly one supersonic root
    f = @(M) area_ratio(M, gamma) - epsilon;
    lo = 1.001; hi = 15.0;
    for i = 1:60   % up to 60 bisection steps -> ~1e-18 precision
        mid = (lo + hi) / 2;
        if f(mid) * f(lo) < 0
            hi = mid;
        else
            lo = mid;
        end
        if (hi - lo) < 1e-10, break; end
    end
    % Use a plain loop instead of for~1 trick for compatibility
    Me = (lo + hi) / 2;
    % Refine with a few more steps
    for k = 1:80
        mid = (lo + hi) / 2;
        if (hi - lo) < 1e-12, break; end
        if f(mid) * f(lo) < 0
            hi = mid;
        else
            lo = mid;
        end
    end
    Me = (lo + hi) / 2;
end


% ── Area ratio A/At as function of Mach ───────────────────────────────
function AR = area_ratio(M, gamma)
    if M <= 0, AR = 1e9; return; end
    t  = 1 + (gamma-1)/2 * M^2;
    AR = (1/M) * (t / ((gamma+1)/2))^((gamma+1)/(2*(gamma-1)));
end


% ── Thrust coefficient — Me passed in, no fsolve ──────────────────────
function Cf = cf_from_me(Pc, Pa, gamma, Me, epsilon)
    % Exit pressure from isentropic relation
    Pe = Pc ./ (1 + (gamma-1)/2 * Me^2).^(gamma/(gamma-1));

    % Momentum thrust coefficient
    Cf_mom = sqrt( 2*gamma^2/(gamma-1) * ...
                   (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
                   (1 - (Pe./Pc).^((gamma-1)/gamma)) );

    % Pressure thrust term
    Cf_pres = (Pe - Pa) .* epsilon ./ Pc;

    Cf = max(Cf_mom + Cf_pres, 0);
end
