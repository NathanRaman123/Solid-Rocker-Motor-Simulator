function grain = BATES(Do, Di, L, N, n_steps)
% BATES  Ballistic Test and Evaluation System grain geometry
%
%   Computes burn area Ab(t), propellant volume Vp(t), and web
%   regression for a multi-segment cylindrical-core BATES grain.
%
%   Inputs:
%     Do      - Outer (case) diameter [m]
%     Di      - Initial core (port) diameter [m]
%     L       - Length of one grain segment [m]
%     N       - Number of grain segments [-]
%     n_steps - Number of time/web steps for discretization [-]
%
%   Outputs:
%     grain   - Struct with fields:
%                 .w        web thickness vector [m]
%                 .Ab       burn area vs. web [m^2]
%                 .Vp       propellant volume vs. web [m^2]
%                 .w_max    maximum web thickness [m]
%                 .Do, .Di, .L, .N  (echoed inputs)
%
%   Physics notes:
%     - Burn progresses radially outward (core expands) AND axially
%       inward from both ends simultaneously.
%     - Web thickness w = (Do - Di) / 2  (radial web)
%     - Each segment burns until core meets OD or ends meet.
%     - Ab = N * (pi * d_core * L_eff  +  2 * pi/4 * (Do^2 - d_core^2))
%       where d_core and L_eff evolve with regression distance x.
%
%   Example:
%     grain = BATES(0.075, 0.025, 0.150, 4, 500);
%     grain_viz(grain);

    % ── Input validation ──────────────────────────────────────────────
    assert(Di < Do, 'Core diameter must be less than outer diameter.');
    assert(Di > 0 && Do > 0 && L > 0, 'Dimensions must be positive.');
    assert(N >= 1 && floor(N) == N, 'N must be a positive integer.');

    % ── Geometry limits ───────────────────────────────────────────────
    w_radial = (Do - Di) / 2;          % Radial web (core to OD) [m]
    w_axial  = L / 2;                  % Axial half-web (burns from both ends) [m]
    w_max    = min(w_radial, w_axial); % Burnout occurs at the smaller

    % ── Regression vector ─────────────────────────────────────────────
    % x = regression distance from initial surface [m]
    x = linspace(0, w_max, n_steps)';

    % ── Instantaneous geometry ────────────────────────────────────────
    d_core  = Di + 2*x;                % Core diameter at regression x [m]
    r_core  = d_core / 2;              % Core radius [m]
    L_eff   = L - 2*x;                 % Effective segment length (both ends burn) [m]
    L_eff   = max(L_eff, 0);           % Clamp — cannot go negative

    % ── Burn areas ────────────────────────────────────────────────────
    % Cylindrical surface (inner bore)
    Ab_cyl  = pi .* d_core .* L_eff;

    % Two annular end faces per segment (outer annulus exposed)
    r_outer = Do / 2;
    Ab_ends = 2 * pi .* (r_outer^2 - r_core.^2);

    % Total for all N segments
    Ab = N .* (Ab_cyl + Ab_ends);

    % Zero out after burnout (both ends met)
    Ab(L_eff <= 0) = 0;

    % ── Propellant volume ─────────────────────────────────────────────
    % Volume of one segment = annular cylinder
    Vp_seg = pi .* (r_outer^2 - r_core.^2) .* L_eff;
    Vp_seg = max(Vp_seg, 0);
    Vp = N .* Vp_seg;

    % ── Burn profile classification ───────────────────────────────────
    % BATES is typically progressive (Ab increases) because the
    % cylindrical surface area grows faster than end area shrinks.
    Ab_mid  = Ab(round(end/2));
    Ab_init = Ab(1);
    if Ab_mid > Ab_init * 1.05
        profile_type = 'Progressive';
    elseif Ab_mid < Ab_init * 0.95
        profile_type = 'Regressive';
    else
        profile_type = 'Neutral';
    end

    % ── Package output ────────────────────────────────────────────────
    grain.type         = 'BATES';
    grain.profile_type = profile_type;
    grain.w            = x;
    grain.Ab           = Ab;
    grain.Vp           = Vp;
    grain.w_max        = w_max;
    grain.Do           = Do;
    grain.Di           = Di;
    grain.L            = L;
    grain.N            = N;
    grain.n_steps      = n_steps;

    fprintf('[BATES] OD=%.1fmm  Core=%.1fmm  L=%.1fmm  N=%d  Web=%.1fmm  Profile=%s\n', ...
        Do*1e3, Di*1e3, L*1e3, N, w_max*1e3, profile_type);
end
